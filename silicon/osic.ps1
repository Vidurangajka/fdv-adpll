<#
.SYNOPSIS
  Run the IIC-OSIC-TOOLS container against silicon/ with sky130A selected.

.DESCRIPTION
  One entry point for every tool the silicon flow uses (ngspice, magic,
  netgen, klayout, xschem), so no command has to remember the details:

    * silicon/ is mounted at /work and is the working directory
    * PDK=sky130A is always set -- the image defaults to IHP sg13g2, and a
      verify run that silently used the wrong PDK is how an earlier version of
      layout/verify.sh got a clean DRC on an empty cell
    * the image tag comes from $env:OSIC_IMAGE, else the pinned default below

.EXAMPLE
  .\osic.ps1 bash layout/verify.sh              # batch: DRC + LVS of ramp_sink
  .\osic.ps1 ngspice -b ngspice/tb_ramp_casc.spice
  .\osic.ps1 -Shell                             # interactive bash in /work
  .\osic.ps1 -Gui                               # desktop at http://localhost:8080
  .\osic.ps1 -Stop                              # stop the -Gui container
  .\osic.ps1 -Check                             # tool versions + PDK sanity
#>
[CmdletBinding(DefaultParameterSetName = "Run")]
param(
    [Parameter(ParameterSetName = "Shell")] [switch]$Shell,
    [Parameter(ParameterSetName = "Gui")]   [switch]$Gui,
    [Parameter(ParameterSetName = "Stop")]  [switch]$Stop,
    [Parameter(ParameterSetName = "Check")] [switch]$Check,
    [Parameter(ParameterSetName = "Run", Position = 0, ValueFromRemainingArguments = $true)]
    [string[]]$Command
)

$ErrorActionPreference = "Stop"
$Image = if ($env:OSIC_IMAGE) { $env:OSIC_IMAGE } else { "hpretl/iic-osic-tools:latest" }
$Work = $PSScriptRoot                   # silicon/
$GuiName = "fdv-osic"
$Common = @("-e", "PDK=sky130A", "-v", "${Work}:/work", "-w", "/work")

docker info --format "{{.ServerVersion}}" *> $null
if ($LASTEXITCODE -ne 0) {
    throw "Docker engine not reachable -- start Docker Desktop first."
}

switch ($PSCmdlet.ParameterSetName) {
    "Gui" {
        # The image's own entry point starts a desktop served over noVNC on
        # port 80 inside the container; publish it on 8080 so it does not
        # collide with anything on the host.
        docker rm -f $GuiName *> $null
        docker run -d --name $GuiName @Common -p 8080:80 $Image | Out-Null
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
        Write-Host "Desktop: http://localhost:8080  (default VNC password: abc123)"
        Write-Host "Stop it with: .\osic.ps1 -Stop"
    }
    "Stop" {
        docker rm -f $GuiName | Out-Null
    }
    "Shell" {
        docker run --rm -it @Common $Image --skip bash
    }
    "Check" {
        docker run --rm @Common $Image --skip bash osic_check.sh
    }
    default {
        if (-not $Command) { Get-Help $PSCommandPath; exit 1 }
        docker run --rm @Common $Image --skip @Command
    }
}
exit $LASTEXITCODE
