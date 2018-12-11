import os
def run():
    dirname = os.path.dirname(__file__)
    dirname = os.path.join(dirname,'launch.sh')
    os.system('bash '+dirname)
def test():
    dirname = os.path.dirname(__file__)
    for filename in os.listdir(dirname):
        if filename.startswith('d0'):
            os.system('jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=60 ' +
                      os.path.join(dirname,filename))
    cleanup()

def cleanup():
    dirname = os.path.dirname(__file__)
    for filename in os.listdir(dirname):
        if filename.endswith('nbconvert.ipynb'):
            os.system('rm ' + os.path.join(dirname,filename))