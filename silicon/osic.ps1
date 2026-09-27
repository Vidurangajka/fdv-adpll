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
# No param() block on purpose.  A declared -Check or -Command would capture any
# argument PowerShell can read as a prefix of it -- `bash -c ...` binds -c to
# them and never reaches the container.  An undeclared script gets every
# argument in $args untouched.
$Mode = if ($args.Count -gt 0 -and $args[0] -in "-Shell", "-Gui", "-Stop", "-Check") {
    $args[0].TrimStart("-")
} else { "Run" }
$Command = $args

# Not "Stop": Windows PowerShell 5.1 would turn every line a tool writes to
# stderr -- ngspice warnings, magic chatter -- into a terminating error.
$Image = if ($env:OSIC_IMAGE) { $env:OSIC_IMAGE } else { "hpretl/iic-osic-tools:latest" }
$Work = $PSScriptRoot                   # silicon/
$GuiName = "fdv-osic"
# A memory cap per container: a runaway simulation is OOM-killed on its own
# instead of starving the Docker VM, which hangs the whole engine (it did,
# twice, on an 8 GB machine).  $env:OSIC_MEM overrides it.
$Mem = if ($env:OSIC_MEM) { $env:OSIC_MEM } else { "2200m" }
$Common = @("-e", "PDK=sky130A", "-v", "${Work}:/work", "-w", "/work",
            "--memory", $Mem, "--memory-swap", $Mem)

docker info --format "{{.ServerVersion}}" *> $null
if ($LASTEXITCODE -ne 0) {
    throw "Docker engine not reachable -- start Docker Desktop first."
}

switch ($Mode) {
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
