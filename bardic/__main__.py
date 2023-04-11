from .cli import bardic_parser


def main():
    args = bardic_parser.parse_args()
    func = args.func
    del args.func
    kwargs = vars(args)
    func(**kwargs)


if __name__ == "__main__":
    main()
