"""Console-script entry points for leafcutter2."""


def leafcutter2_cli():
    """Entry point for the ``leafcutter2`` command."""
    from leafcutter2.pipeline import main_cli
    main_cli()


def make_clusters_cli():
    """Entry point for the ``leafcutter2-make-clusters`` command."""
    from leafcutter2.clustering import main_cli
    main_cli()


def transcript_tools_cli():
    """Entry point for the ``leafcutter2-transcript-tools`` command."""
    from leafcutter2.transcript_tools import main
    main()


def star2junc_cli():
    """Entry point for the ``leafcutter2-star2junc`` command."""
    from leafcutter2.star_utils import main_cli
    main_cli()
