from asyncio import Event

import ray
from tqdm import tqdm


@ray.remote
class ProgressBarActor:

    def __init__(self):
        self.counter = 0
        self.delta = 0
        self.event = Event()

    def update(self, num_items_completed):
        """Updates the ProgressBar with the incremental
        number of items that were just completed.
        """
        self.counter += num_items_completed
        self.delta += num_items_completed
        self.event.set()

    async def wait_for_update(self):
        """Blocking call.

        Waits until somebody calls `update`, then returns a tuple of
        the number of updates since the last call to
        `wait_for_update`, and the total number of completed items.
        """
        await self.event.wait()
        self.event.clear()
        saved_delta = self.delta
        self.delta = 0
        return saved_delta, self.counter

    def get_counter(self):
        """
        Returns the total number of complete items.
        """
        return self.counter


class ProgressBar:

    def __init__(self, total, description=""):
        # Ray actors don't seem to play nice with mypy, generating
        # a spurious warning for the following line,
        # which we need to suppress. The code is fine.
        self.progress_actor = ProgressBarActor.remote()  # type: ignore
        self.total = total
        self.description = description

    @property
    def actor(self):
        """Returns a reference to the remote `ProgressBarActor`.

        When you complete tasks, call `update` on the actor.
        """
        return self.progress_actor

    def print_until_done(self) -> None:
        """Blocking call.

        Do this after starting a series of remote Ray tasks, to which you've
        passed the actor handle. Each of them calls `update` on the actor.
        When the progress meter reaches 100%, this method returns.
        """
        pbar = tqdm(desc=self.description, total=self.total)
        while True:
            delta, counter = ray.get(self.actor.wait_for_update.remote())
            pbar.update(delta)
            if counter >= self.total:
                pbar.close()
                return


def process_by_chunks(iterable, remote_func, chunksize, *args, **kwargs):
    """Parallelizes function execution with ray in chunks of iterable.
    
    Args:
        iterable (iterable): an iterable with __len__ attribute.
            The elements are passed one by one into the function
            as the first argument.
        remote_func (ray.remote_function.RemoteFunction): a ray remote function
            to be executed.
        chunksize (int): a number of simultanious task to schedule.
        args: position arguments to be passed to the function.
        kwargs: keyword arguments to be passed to the function.
    
    Returns:
        list: function execution results.
    """
    chunk_indices = list(range(0, len(iterable), chunksize))
    if chunk_indices[-1] < len(iterable):
        chunk_indices.append(len(iterable))
    results = list()
    for i, (start, end) in tqdm(enumerate(zip(chunk_indices[:-1], chunk_indices[1:])), desc='Chunks'):
        chunk = iterable[start:end]
        pb = ProgressBar(len(chunk))
        tasks_pre_launch = [remote_func.remote(item, *args, pb.actor, **kwargs) for item in chunk]
        pb.print_until_done()
        chunk_results = ray.get(tasks_pre_launch)
        results += chunk_results
    return results
