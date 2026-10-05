@echo off
setlocal
powershell.exe -NoLogo -NoProfile -NonInteractive -ExecutionPolicy Bypass -File "%~dp0launch.ps1" %*
set "launchExitCode=%ERRORLEVEL%"
exit /b %launchExitCode%