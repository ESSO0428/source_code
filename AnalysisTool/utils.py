"""
    Python utils for the data analysis in `AnalysisTool`
"""

from typing import Protocol, TypeVar

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


_T_contra = TypeVar("_T_contra", contravariant=True)
# stable
class SupportsWrite(Protocol[_T_contra]):
    def write(self, __s: _T_contra) -> object: ...


def tee(data, file: SupportsWrite[str] | None, end='\n'):
    """
    Prints and writes to a file at the same time.
    :param data: The data to be printed and written to the file.
    :param file: The file object to which the data should be written.
    :param end: The character(s) to append at the end. Default is newline.
    """
    print(data, end=end)  # prints to console
    if file:
        print(data, end=end, file=file)  # writes to file

