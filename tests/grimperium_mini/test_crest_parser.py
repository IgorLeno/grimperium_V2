from grimperium_mini.crest import split_crest_conformers


def test_split_crest_conformers(tmp_path):
    ensemble = tmp_path / "crest_conformers.xyz"
    ensemble.write_text(
        "2\nconf1\nH 0 0 0\nH 0 0 1\n\n2\nconf2\nH 0 0 0\nH 0 1 0\n",
        encoding="utf-8",
    )
    files = split_crest_conformers(ensemble, tmp_path / "split")
    assert [p.name for p in files] == ["conformer_0001.xyz", "conformer_0002.xyz"]
    assert files[1].read_text(encoding="utf-8").splitlines()[1] == "conf2"
