"""
    Python utils for the data analysis in `AnalysisTool`
"""


def iter_status(iterable, start=0):
    """Pass through all values from the given iterable, starting with index 'start' (default 0).
    Each item is a tuple containing:
    1. The index of the item.
    2. A boolean which is True if there are more items to come, and False if it is the last item.
    3. The item itself.
    ---

    for example:
    >>> iterable = ['A', 'B', 'C']
    >>> for i, has_next, value in iter_status(iterable):
            print(i, has_next, value)
    0 True A
    1 True B
    2 False C 
    """
    it = iter(iterable)
    index = start
    last_item = next(it)
    for item in it:
        yield index, True, last_item
        index += 1
        last_item = item
    yield index, False, last_item

