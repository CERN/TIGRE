import os
def test():
    """
    Function only outputs to terminal but gives good coverage of
    the functions implemented and works well as a unit test for
    imports.

    Usage
    --------
    import tigre
    tigre.demos.test()

    :return: None
    """
    dirname = os.path.dirname(__file__)
    for filename in os.listdir(dirname):
        if filename.startswith('d0'):
            """
            cmd = ['jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=60 ' +
                      os.path.join(dirname,filename)]
            print(cmd)
            output=subprocess.Popen(cmd, stdout=subprocess.PIPE,shell=True).communicate()[0]
            nb = nbformat.read(output, nbformat.current_nbformat)
            """
            os.system('jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=60 ' +
                      os.path.join(dirname,filename))

    cleanup()

def cleanup():
    dirname = os.path.dirname(__file__)
    for filename in os.listdir(dirname):
        if filename.endswith('nbconvert.ipynb'):
            os.system('rm ' + os.path.join(dirname,filename))