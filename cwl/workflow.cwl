cwlVersion: v1.0
class: Workflow
inputs:
  url: string

outputs:
  ms:
    type: Directory
    outputSource: h5toms/ms


steps:
  curl:
    run: curl.cwl
    in:
      url: url
    out: [downloaded]

  h5toms:
    run: h5toms.cwl
    in:
      h5: curl/downloaded
    out: [ms]