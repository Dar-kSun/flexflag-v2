"""Command-line entry point. `flexflag check` lands in M7."""

import typer

from flexflag import __version__

app = typer.Typer(help="Flag proteins likely to change shape between apo and holo states.")


@app.callback()
def main() -> None:
    """Keep subcommands explicit even while there is only one."""


@app.command()
def version() -> None:
    """Print the installed version."""
    typer.echo(__version__)


if __name__ == "__main__":
    app()
