import numpy as np
from typing import List, Tuple, Dict


def einsum_with_names(
    terms: List[Tuple[np.ndarray, Tuple[str, ...]]],
    output_labels: Tuple[str, ...],
) -> np.ndarray:
    """
    terms: list of (array, labels_tuple) where labels_tuple are composite label strings
           e.g. ("i(m,q)","i(m+1,0)","a(m+1,1)","p(3)")
    output_labels: tuple of composite labels for the result (order matters)
    Returns: np.einsum result using the integer-index API so label names can be arbitrary strings.
    """
    # collect all unique labels and assign integer ids
    label_to_int: Dict[str, int] = {}
    next_int = 0

    def get_int(lbl):
        nonlocal next_int
        if lbl not in label_to_int:
            label_to_int[lbl] = next_int
            next_int += 1
        return label_to_int[lbl]

    arrays = []
    index_lists = []

    # convert each term's label tuple into a list of integers
    for arr, labels in terms:
        idxs = [get_int(lbl) for lbl in labels]
        arrays.append(arr)
        index_lists.append(idxs)

    # ensure output labels are in the mapping
    out_idxs = [get_int(lbl) for lbl in output_labels]

    # build the einsum call: np.einsum(arr0, idx0, arr1, idx1, ..., out_idx_list)
    einsum_args = []
    for arr, idxs in zip(arrays, index_lists):
        einsum_args.append(arr)
        einsum_args.append(idxs)
    einsum_args.append(out_idxs)

    return np.einsum(*einsum_args)


def gen_state(comb: List[int], nb_qudits: int, qudit_len: int):
    state = {"qudits": [], "roots": []}

    for i, label in enumerate(comb):
        if i < nb_qudits * qudit_len:
            if i % qudit_len:
                state["qudits"][-1].append(label)
            else:
                state["qudits"].append([label])
        else:
            state["roots"].append(label)

    return state
