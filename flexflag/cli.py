"""Command-line entry point: `flexflag check <UniProt> --residues ... | --from-pdb ...`."""

import typer

from flexflag import __version__
from flexflag.check import check, format_report, parse_residues, pocket_from_pdb

app = typer.Typer(help="Flag protein pockets likely to change shape between apo and holo states.")


@app.callback()
def main() -> None:
    """Keep subcommands explicit."""


@app.command("check")
def check_cmd(
    uniprot: str = typer.Argument(..., help="UniProt accession, e.g. P69441"),
    residues: str = typer.Option(
        None, "--residues", "-r", help="Pocket residues, UniProt numbering: '13,31,35-38'"
    ),
    from_pdb: str = typer.Option(
        None, "--from-pdb", help="PDB ID or local mmCIF/PDB file with a bound ligand"
    ),
) -> None:
    """Risk that the pocket changes shape between empty and bound, from AlphaFold DB."""
    if residues and from_pdb:
        raise typer.BadParameter("give --residues or --from-pdb, not both")
    uniprot = uniprot.strip().upper()
    pocket, source = None, None
    try:
        if residues:
            pocket = parse_residues(residues)
        elif from_pdb:
            pocket, source = pocket_from_pdb(uniprot, from_pdb)
        report = check(uniprot, pocket, source)
    except (ValueError, LookupError) as e:
        typer.echo(f"error: {e}", err=True)
        raise typer.Exit(1) from e
    typer.echo(format_report(report))


@app.command()
def version() -> None:
    """Print the installed version."""
    typer.echo(__version__)


if __name__ == "__main__":
    app()
