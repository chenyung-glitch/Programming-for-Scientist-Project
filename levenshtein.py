"""Utilities for comparing sequences by edit distance."""


def levenshtein_distance(first: str, second: str) -> int:
    """Return the minimum edits needed to turn ``first`` into ``second``.

    An edit is one insertion, deletion, or substitution. The implementation
    keeps only the previous row of the dynamic-programming table, so its
    memory use grows with the length of the shorter input.
    """
    if len(first) < len(second):
        first, second = second, first

    previous_row = list(range(len(second) + 1))

    for first_index, first_character in enumerate(first, start=1):
        current_row = [first_index]
        for second_index, second_character in enumerate(second, start=1):
            insertion = current_row[second_index - 1] + 1
            deletion = previous_row[second_index] + 1
            substitution = previous_row[second_index - 1] + (
                first_character != second_character
            )
            current_row.append(min(insertion, deletion, substitution))
        previous_row = current_row

    return previous_row[-1]
