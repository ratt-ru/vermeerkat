cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: vermeerkat/curl
baseCommand: [curl, -O, -L, -C, -]
inputs:
  url:
    type: string
    inputBinding:
      position: 1
outputs:
  h5:
    type: File
    outputBinding:
      glob: "*"
