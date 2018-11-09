#!/usr/bin/env cwl-runner
#
# Example validate submission file
#
cwlVersion: v1.0
class: CommandLineTool
baseCommand: python

inputs:
  - id: inputfile
    type: File
  - id: goldstandard
    type: File

arguments:
  - valueFrom: validate.py
  - valueFrom: $(inputs.inputfile)
    prefix: -s
  - valueFrom: results.json
    prefix: -r
  - valueFrom: $(inputs.goldstandard.path)
    prefix: -g

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: validate.py
        entry: |
          #!/usr/bin/env python
          import synapseclient
          import argparse
          import os
          import json
          import pandas as pd
          parser = argparse.ArgumentParser()
          parser.add_argument("-s", "--submission_file", required=True, help="Submission File")
          parser.add_argument("-r", "--results", required=True, help="validation results")
          parser.add_argument("-g", "--goldstandard", required=True, help="Goldstandard for scoring")

          args = parser.parse_args()

          def _validate_func_helper(filePath, goldDf, predOrConf, column="proteinID", varianceCheck=False, scoring_sc1=True):
            fileName = os.path.basename(filePath)
            if scoring_sc1:
              assert os.path.isfile(filePath), "%s file must be named %s, and your model must generate %s_{1..100}.tsv" % (predOrConf, fileName, predOrConf)
            else:
              assert os.path.isfile(filePath), "%s file must be named %s, and your model must generate %s.tsv" % (predOrConf, fileName, predOrConf)

            assert os.stat(filePath).st_size > 0, "%s: Can't be an empty file" % fileName
            try :
              fileDf = pd.read_csv(filePath, sep="\t",header=None)
              assert fileDf[0][0] == column, "Please do not write out the row names in your prediction file, and %s must be the first column header" % column
            except pd.errors.ParserError as e:
              raise AssertionError("Please do not write out the row names in your prediction file.")

            fileDf = pd.read_csv(filePath, sep="\t")

            assert fileDf.get(column) is not None, "%s: Must contain %s column" % (fileName,column)
            fileDf.index = fileDf[column]
            del fileDf[column]
            assert all(~fileDf.index.duplicated()), "%s: There cannot be any duplicated %s" % (fileName,column)
            assert all(~fileDf.columns.duplicated()), "%s: There cannot be any duplicated sample ids" % fileName
            assert all(goldDf.index.isin(fileDf.index)), "%s: All %s in the goldstandard must also be in your file. You are missing: %s" % (fileName, column, ",".join(set(goldDf.index[~goldDf.index.isin(fileDf.index)].map(str))))
            assert all(goldDf.columns.isin(fileDf.columns)), "%s: All sample Ids in the goldstandard must also be in your file. You are missing: %s" % (fileName, ",".join(goldDf.columns[~goldDf.columns.isin(fileDf.columns)]))
            assert sum(fileDf.isnull().sum()) == 0, "%s: There can't be any null values" % fileName
            assert all(fileDf.applymap(lambda x: isinstance(x, (float,int))).all()), "%s: Must submit float or int values" % fileName

          invalid_reasons = []
          prediction_file_status = "VALIDATED"
          validation_message = None
          golddf = pd.read_csv(args.goldstandard, sep="\t",index_col=0)

          try:
            _validate_func_helper(args.submission_file, golddf, "predictions", column="proteinID", varianceCheck=True, scoring_sc1=False)
          except AssertionError as ex1:
            validation_message = str(ex1)

          if validation_message is not None:
            invalid_reasons.append(validation_message)
            prediction_file_status = "INVALID"
          result = {'prediction_file_errors':"\n".join(invalid_reasons),'prediction_file_status':prediction_file_status}
          with open(args.results, 'w') as o:
            o.write(json.dumps(result))
     
outputs:

  - id: results
    type: File
    outputBinding:
      glob: results.json   

  - id: status
    type: string
    outputBinding:
      glob: results.json
      loadContents: true
      outputEval: $(JSON.parse(self[0].contents)['prediction_file_status'])

  - id: invalid_reasons
    type: string
    outputBinding:
      glob: results.json
      loadContents: true
      outputEval: $(JSON.parse(self[0].contents)['prediction_file_errors'])
