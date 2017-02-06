# from bitarray import bitarray
import numpy as np

__all__ = ['BinarySubspace']
def xnor(a,b):
    return (a&b)^(~a&~b)


def xor(a,b):
      return (a|b)&~(a&b)


class BinarySubspace(object):
    """Set-like class for bitarray objects to generate a closed subspace."""
    def __init__(self, *data):
        self.order = 0
        self._items = []
        self.generators = []
        for val in data:
            # if not isinstance(val, bitarray):
            if not isinstance(val, np.ndarray):
                raise ValueError('This class works for numpy arrays only!')
                # raise ValueError('This class works for bitarrays only!')
            self.add(val)
    
    def __contains__(self, it):
        for _el in self._items:
            # if all(xnor(_el, it)):
            if np.array_equal(_el, it):
                return True
        return False

    def __iter__(self):
        for item in self._items:
            yield item

    def _generate(self, obj):
        for item in self._items:
            # new = item^obj
            new = xor(item, obj)
            if new in self:
                continue
            else:
                self.order +=1
                self._items.append(new)
                self._generate(new)
        return

    def __eq__(self, other):
        return all([_el in other for _el in self._items])

    def add(self, obj):
        for _el in self._items:
            if all(xnor(obj, _el)):
                return self
        self.order +=1
        self.generators.append(obj)
        self._items.append(obj)
        self._generate(obj)
        return self