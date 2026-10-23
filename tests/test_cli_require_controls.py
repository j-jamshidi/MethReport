"""--require-controls turns a missing control file into an error."""

from typer.testing import CliRunner

from methreport.cli import app

runner = CliRunner()


def _args(tmp_path, *extra):
    bed = tmp_path / "unphased.bed"
    bed.write_text("")
    return ["run", "--bed", str(bed), "--ref", "hg38", "--out", str(tmp_path / "out"), *extra]


def test_fails_when_no_controls_are_available(tmp_path):
    result = runner.invoke(app, _args(tmp_path, "--require-controls"))
    assert result.exit_code == 1
    assert "no control data was loaded" in result.output.replace("\n", " ")


def test_fails_on_a_missing_controls_file(tmp_path):
    result = runner.invoke(
        app, _args(tmp_path, "--require-controls", "--controls", str(tmp_path / "nope.tsv")))
    assert result.exit_code == 1
    assert "could not be loaded" in result.output.replace("\n", " ")
