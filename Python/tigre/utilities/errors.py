class Error(Exception):
    """Base class for exceptions in this module"""

    def __init__(self):
        Exception.__init__(self)


class TigreCudaCallError(Error):
    """Error list indexed against values in errors.hpp."""

    def __init__(self, txt, error_index):
        Error.__init__(self)
        self.error_list = [
            "GPU: Cuda error with exit code 1.",
            "There are no available device(s) that support CUDA",
            "GPU: One (or more) of your GPUs is being heavily used by another program (possibly graphics-based).\n"  # noqa: E501
            "Free the GPU to run TIGRE",
            "GPU: Assertion Failed. Logic behind spliting flawed! Please tell: ander.biguri@gmail.com",  # noqa: E501
            "GPU: Assertion Failed. Memory needed calculation broken! Please tell: ander.biguri@gmail.com",  # noqa: E501
        ]
        self.error_index = error_index
        self.txt = txt

    def __str__(self):
        return "\n".join([self.txt, self.error_list[self.error_index - 1]])
