cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: vermeerkat/h5toms
baseCommand: [h5toms.py]
inputs:
  h5:
    type: file
    inputBinding:
      position: 1
outputs:
  ms:
    type: Directory
    outputBinding:
      glob: "*"
