cwlVersion: v1.0
class: Workflow
inputs:
  url: string

outputs:
  ms:
    type: Directory
    outputSource: ms

steps:
  curl:
    run: curl.cwl
    in:
      url: url
    out: [h5]

  compile:
    run: h5toms.cwl
    in:
      h5: h5
    out: [ms]