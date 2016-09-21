__all__ = ['format_for_output', 'AnalysisResult']

class AnalysisResult(object):
    def __init__(self, n_qubits, target, *args, **kwargs):
        self.n = n_qubits
        self.target = target
        self.fname = kwargs.pop('fname', None)
        self.ostring = kwargs.pop('ostring', '')
        self.data = args
    def __repr__(self):
        return self.ostring.format(n_qubits=self.n_qubits, 
                                   target=self.target,
                                   *self.data)
    def write(self):
        if self.fname is None:
            print(str(self))
        else:
            with open(self.fname, 'a') as f:
                f.write(str(self))

