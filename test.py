from levenshtein import levenshtein_distance


def test_levenshtein_distance():
    assert levenshtein_distance("kitten", "sitting") == 3
    assert levenshtein_distance("GATTACA", "GACTATA") == 2


def test_levenshtein_distance_handles_empty_strings():
    assert levenshtein_distance("", "") == 0
    assert levenshtein_distance("abc", "") == 3
    assert levenshtein_distance("", "abc") == 3


def test_levenshtein_distance_handles_identical_strings():
    assert levenshtein_distance("genome", "genome") == 0
