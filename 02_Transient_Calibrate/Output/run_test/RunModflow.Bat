Copy /Y /A "Transient_cali.mfsim.nam" "mfsim.nam"
call "D:\Program Files\USGS\ModelMuse5\bin\ModelMonitor.exe" -m "D:\Program Files\USGS\mf6.6.2_win64\bin\mf6.exe" -n mfsim.nam -mv 6
Start Notepad.exe Transient_cali.lst
Start Notepad.exe mfsim.lst
pause
