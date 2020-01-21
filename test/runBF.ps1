# Change the 2 to whatever 2 <= n <= 30 you want to test
For ($j = 0; $j -lt 10; $j++) {
    $fileAsStr = ".\Inputs\input_2_$j"
    $fileAsDir = Get-Item $fileAsStr
    # echo $fileAsDir
    Measure-Command { Get-Content $fileAsDir | .\SetCoverBF\SetCoverBF.exe } | Select-Object -Property TotalSeconds
}