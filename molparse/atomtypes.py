from ase.data import chemical_symbols

chemical_symbols = chemical_symbols[1:]
chemical_symbols = sorted(chemical_symbols, key=lambda x: len(x), reverse=True)


def formula_to_atomtype_dict(formula):
    """

    :param formula:

    """

    import re

    remaining = formula.replace("+", "").replace("-", "")

    atomtype_dict = {}

    while remaining:
        for symbol in chemical_symbols:
            if remaining.startswith(symbol):
                remaining = remaining.removeprefix(symbol)

                count_str = re.findall(r"^([0-9]*)", remaining)[0]

                if count_str:
                    remaining = remaining.removeprefix(count_str)
                    count = int(count_str)
                else:
                    count = 1

                atomtype_dict[symbol] = count

                break
        else:
            raise Exception(
                f"Could not match any chemical symbol to the start of {remaining} {formula=}"
            )

    return atomtype_dict


def combine_atomtype_dicts(atomtype_dicts):
    """

    :param atomtype_dicts:

    """

    merged = {}
    for atomtype_dict in atomtype_dicts:
        for symbol, count in atomtype_dict.items():
            if symbol not in merged:
                merged[symbol] = 0

            merged[symbol] += count

    return merged


def atomtype_dict_to_formula(atomtype_dict):
    """

    :param atomtype_dict:

    """

    symbols = []

    """For organic compounds, the order is carbon, hydrogen,
	then all other elements in alphabetical order of their chemical symbols."""

    key = "C"
    if key in atomtype_dict:
        value = atomtype_dict[key]
        if value > 1:
            symbols.append(f"C{value}")
        else:
            symbols.append("C")

    key = "H"
    if key in atomtype_dict:
        value = atomtype_dict[key]
        if value > 1:
            symbols.append(f"{key}{value}")
        else:
            symbols.append(key)

    keys = sorted(list(atomtype_dict.keys()))

    for key in keys:
        value = atomtype_dict[key]
        if key in ["C", "H"]:
            continue
        if value > 1:
            symbols.append(f"{key}{value}")
        else:
            symbols.append(key)

    return "".join(symbols)


def subtract_atomtype_dict(d1, d2, ignore_hydrogen=False):

    all_keys = set(d1.keys()) | set(d2.keys())

    d = {}

    for key in all_keys:

        if ignore_hydrogen and key == "H":
            continue

        n1 = d1[key] if key in d1 else 0
        n2 = d2[key] if key in d2 else 0

        n = n1 - n2

        if n > 0:
            d[key] = n

    return d
