from typing import Optional


def adjust_chunksize(iterable_size: int, n_cores: int, chunksize: Optional[int] = None, limit: int = 1000):
    """Adjusts a user-provided chunk size to the data and a number of workers.

    If a size of an iterable is larger than specified value, the data should be
    split into chunks larger than 1 to increase the performance. A data-tailored
    chunk size is calculated as a ratio between the size of the iterable and a
    number of cores.

    User can provide their own chunk size. In that case, adjusted chunk size
    is a minimum of the user-provided and data-tailored chunk sizes.

    Parameters
    ----------
    iterable_size
        The size of the iterable that will be passed to a ProcessPoolExecutor instance.
    n_cores
        Number of cores used (max_workers of the executor instance).
    chunksize
        The user-provided chunk size.
    limit
        A maximal size of the iterable to set the chunk size to 1.

    Returns
    -------
    int
        The adjusted chunk size.
    """
    if iterable_size < limit:
        return 1
    data_chunks = max(1, iterable_size // n_cores)
    if chunksize is None:
        return data_chunks
    return min(data_chunks, chunksize)
