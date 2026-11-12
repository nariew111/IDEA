from idea.cli import main
from idea_gp.cli import main as gp_main
import idea_gp.cli as gp_cli
import pytest


def test_top_level_help(capsys):
    try:
        main(["--help"])
    except SystemExit as exc:
        assert exc.code == 0
    out = capsys.readouterr().out
    assert "idea-design" in out
    assert "generate-features" in out
    assert "select" in out


def test_idea_gp_status(capsys):
    gp_main(["status"])
    out = capsys.readouterr().out
    assert "HAGSSuite" in out
    assert "Runner attached: True" in out


@pytest.mark.parametrize("separator", [[], ["--"]])
def test_idea_gp_forwards_dataset_and_output_options(monkeypatch, separator):
    received = []
    monkeypatch.setattr(gp_cli, "_run_hags", received.append)
    runner_args = ["--bfile", "genotypes", "--pheno", "traits.tsv", "--out", "results/gp"]
    gp_main(["run", *separator, *runner_args])
    assert received == [runner_args]


def test_idea_gp_help_reaches_the_runner(monkeypatch):
    received = []
    monkeypatch.setattr(gp_cli, "_run_hags", received.append)
    gp_main(["run", "--help"])
    assert received == [["--help"]]
