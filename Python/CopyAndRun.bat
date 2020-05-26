xcopy build\lib.win-amd64-3.7\_Atb.cp37-win_amd64.pyd         .\ /y
xcopy build\lib.win-amd64-3.7\_AwminTV.cp37-win_amd64.pyd     .\ /y
xcopy build\lib.win-amd64-3.7\_Ax.cp37-win_amd64.pyd          .\ /y
xcopy build\lib.win-amd64-3.7\_gpuUtils.cp37-win_amd64.pyd    .\ /y
xcopy build\lib.win-amd64-3.7\_minTV.cp37-win_amd64.pyd       .\ /y
xcopy build\lib.win-amd64-3.7\_tvdenoising.cp37-win_amd64.pyd .\ /y
rem xcopy build\lib.win-amd64-3.7\_gpuUtils.cp37-win_amd64.pyd    .\ /y

Python tests\list_gpu.py
REM pause
python example.py

