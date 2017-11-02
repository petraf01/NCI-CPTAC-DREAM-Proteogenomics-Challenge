# Use rpy2 if you have R scoring functions
import rpy2.robjects as robjects
import os
import zipfile
import pandas as pd
filePath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'scoring_functions.R')
robjects.r("source('%s')" % filePath)
corr_by_row = robjects.r('correlation_by_row')
nrmse_by_row = robjects.r('NRMSE_by_row')

score_cor = robjects.r('score.cor')
score_nrmse = robjects.r('score.nrmsd')

corr_by_row_sc3 = robjects.r('correlation_by_row_ALL_OBSERVED')
nrmse_by_row_sc3 = robjects.r('NRMSE_by_row_ALL_OBSERVED')

##-----------------------------------------------------------------------------
##
## challenge specific code and configuration
##
##-----------------------------------------------------------------------------


## A Synapse project will hold the assetts for your challenge. Put its
## synapse ID here, for example
## CHALLENGE_SYN_ID = "syn1234567"
CHALLENGE_SYN_ID = "syn8228300"

## Name of your challenge, defaults to the name of the challenge's project
CHALLENGE_NAME = "NCI-CPTAC Proteogenomics Challenge"

## Synapse user IDs of the challenge admins who will be notified by email
## about errors in the scoring script
ADMIN_USER_IDS = [3324230,3360851]


## Each question in your challenge should have an evaluation queue through
## which participants can submit their predictions or models. The queues
## should specify the challenge project as their content source. Queues
## can be created like so:
##   evaluation = syn.store(Evaluation(
##     name="My Challenge Q1",
##     description="Predict all the things!",
##     contentSource="syn1234567"))
## ...and found like this:
##   evaluations = list(syn.getEvaluationByContentSource('syn3375314'))
## Configuring them here as a list will save a round-trip to the server
## every time the script starts and you can link the challenge queues to
## the correct scoring/validation functions.  Predictions will be validated and 

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
    # if varianceCheck:
    #     assert all(fileDf.apply(pd.np.var,axis=1) != 0), "%s: No rows can have variance of 0" % fileName


def validate_func1(dirName, goldstandardDir, column):
    for num in range(1,101):
        goldstandard_path = os.path.join(goldstandardDir, "data_test_true_%s.txt" % num)
        goldDf = pd.read_csv(goldstandard_path, sep="\t",index_col=0)   
        prediction_path = os.path.join(dirName,'predictions_%d.tsv' % num)
        confidence_path = os.path.join(dirName,'confidence_%d.tsv' % num)
        _validate_func_helper(prediction_path, goldDf, "predictions")
        _validate_func_helper(confidence_path, goldDf, "confidence")
    return(True,"Passed Validation")

def validate_func2_3(dirName, goldstandard_path, column):
    if column = "proteinID"
        breast_goldDf = pd.read_csv(os.path.join(goldstandard_path,"prospective_breast_proteome_sort_common_gene_10005.txt"), sep="\t",index_col=0)
        ova_goldDf = pd.read_csv(os.path.join(goldstandard_path,"prospective_ova_proteome_sort_common_gene_7061.txt"), sep="\t",index_col=0)
    else:
        breast_goldDf = pd.read_csv(os.path.join(goldstandard_path,"prospective_breast_phospho_sort_common_gene_31981.txt"), sep="\t",index_col=0)
        ova_goldDf = pd.read_csv(os.path.join(goldstandard_path,"prospective_ova_phospho_sort_common_gene_10057.txt"), sep="\t",index_col=0)
       
    breast_prediction_path = os.path.join(dirName,'breast_predictions.tsv')
    breast_confidence_path = os.path.join(dirName,'breast_confidence.tsv')
    ovarian_prediction_path = os.path.join(dirName, 'ovarian_predictions.tsv')
    ovarian_confidence_path = os.path.join(dirName,'ovarian_confidence.tsv')

    if os.path.exists(breast_prediction_path):
        _validate_func_helper(breast_prediction_path, goldDf, "predictions", column=column, varianceCheck=True, scoring_sc1=False)
        _validate_func_helper(breast_confidence_path, goldDf, "confidence", column=column, scoring_sc1=False)
        _validate_func_helper(ovarian_prediction_path, goldDf, "predictions", column=column, varianceCheck=True, scoring_sc1=False)
        _validate_func_helper(ovarian_confidence_path, goldDf, "confidence", column=column, scoring_sc1=False)
    else:
        goldDf = pd.read_csv(goldstandard_path, sep="\t",index_col=0)   
        breast_prediction_path = os.path.join(dirName,'predictions.tsv')
        _validate_func_helper(breast_prediction_path, goldDf, "predictions", column=column, varianceCheck=True, scoring_sc1=False)
        _validate_func_helper(breast_confidence_path, goldDf, "confidence", column=column, scoring_sc1=False)

    return(True,"Passed Validation")

def score1(dirName, goldstandardDir):
    ##Read in submission (submission.filePath)
    ##Score against goldstandard
    #goldstandard_path = os.path.join(goldstandardDir, "data_test_true.txt")
    sc1_nrmsd_scores = []
    sc1_corr_scores = []
    for num in range(1,101):
        goldstandard_path = os.path.join(goldstandardDir, "data_test_true_%s.txt" % num)
        prediction_path = os.path.join(dirName,'predictions_%d.tsv' % num)
        observed_path = os.path.join(goldstandardDir, "data_test_obs_%d.txt" % num)
        sc1_corr_scores.append(score_cor(prediction_path, observed_path, goldstandard_path)[0])
        sc1_nrmsd_scores.append(score_nrmse(prediction_path, observed_path, goldstandard_path)[0])

    return(pd.np.mean(sc1_corr_scores), pd.np.mean(sc1_nrmsd_scores))

def score2(dirName, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Score against goldstandard

    prediction_path = os.path.join(dirName,'predictions.tsv')
    if os.path.exists(prediction_path):
        corr = corr_by_row(prediction_path, goldstandard_path)[0]
        rmse = nrmse_by_row(prediction_path, goldstandard_path)[0]
        preds = dict(corr=round(corr,4), rmse=round(rmse,4))
    else:
        breast_goldPath = os.path.join(goldstandard_path,"prospective_breast_proteome_sort_common_gene_10005.txt")
        ova_goldPath = os.path.join(goldstandard_path,"prospective_ova_proteome_sort_common_gene_7061.txt")
        breast_prediction_path = os.path.join(dirName,'breast_predictions.tsv')
        breast_confidence_path = os.path.join(dirName,'breast_confidence.tsv')
        ovarian_prediction_path = os.path.join(dirName, 'ovarian_predictions.tsv')
        ovarian_confidence_path = os.path.join(dirName,'ovarian_confidence.tsv')
        breast_corr = corr_by_row(breast_prediction_path, breast_goldPath)[0]
        breast_rmse = nrmse_by_row(breast_prediction_path, breast_goldPath)[0] 
        ova_corr = corr_by_row(ovarian_prediction_path, ova_goldPath)[0]
        ova_rmse = nrmse_by_row(ovarian_prediction_path, ova_goldPath)[0] 
        preds = dict(breast_corr=round(breast_corr,4), breast_rmse=round(breast_rmse,4), ova_corr=round(ova_corr,4), ova_rmse=round(ova_rmse,4))
    return(preds)

def score3(dirName, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Score against goldstandard
    prediction_path = os.path.join(dirName,'predictions.tsv')
    if os.path.exists(prediction_path):
        corr = corr_by_row(prediction_path, goldstandard_path)[0]
        rmse = nrmse_by_row(prediction_path, goldstandard_path)[0]
        preds = dict(corr=round(corr,4), rmse=round(rmse,4))
    else:
        breast_goldPath = os.path.join(goldstandard_path,"prospective_breast_phospho_sort_common_gene_31981.txt")
        ova_goldPath = os.path.join(goldstandard_path,"prospective_ova_phospho_sort_common_gene_10057.txt")
        breast_prediction_path = os.path.join(dirName,'breast_predictions.tsv')
        breast_confidence_path = os.path.join(dirName,'breast_confidence.tsv')
        ovarian_prediction_path = os.path.join(dirName, 'ovarian_predictions.tsv')
        ovarian_confidence_path = os.path.join(dirName,'ovarian_confidence.tsv')
        breast_corr = corr_by_row(breast_prediction_path, breast_goldPath)[0]
        breast_rmse = nrmse_by_row(breast_prediction_path, breast_goldPath)[0] 
        ova_corr = corr_by_row(ovarian_prediction_path, ova_goldPath)[0]
        ova_rmse = nrmse_by_row(ovarian_prediction_path, ova_goldPath)[0] 
        preds = dict(breast_corr=round(breast_corr,4), breast_rmse=round(breast_rmse,4), ova_corr=round(ova_corr,4), ova_rmse=round(ova_rmse,4))
    return(preds)
#3 weeks
#quota =  
# {u'firstRoundStart': u'2017-07-14T00:00:00.000Z',
#   u'numberOfRounds': 1,
#   u'roundDurationMillis': 3369599000,
#   u'submissionLimit': 3}
#1503385199000 - 1500015600000
evaluation_queues = [

# Proteogenomics Subchallenge 1 (8720143)
# Proteogenomics Subchallenge 2 (8720145)
# Proteogenomics Subchallenge 3 (8720149)
    {
        'id':8720143,
        'scoring_func':score1,
        'validation_func':validate_func1,
        'column':'proteinID',
        #THIS IS AN EXCEPTION.  NEEDS TO PROVIDE GOLDSTANDARD DIRECTORY
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'goldstandard')
    },
    {
        'id':8720145,
        'scoring_func':score2,
        'validation_func':validate_func2_3,
        'column':'proteinID',
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'goldstandard/prospective_ova_pro_gold.txt')
    },
    {
        'id':8720149,
        'scoring_func':score3,
        'validation_func':validate_func2_3,
        'column':'phosphoID',
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'goldstandard/prospective_ova_phospho_gold.txt')
    },
# Proteogenomics Subchallenge 1 Express Lane (9604716)
# Proteogenomics Subchallenge 2 Express Lane (9604717)
# Proteogenomics Subchallenge 3 Express Lane (9604718)
    {
        'id':9604716,
        'scoring_func':None,
        'validation_func':validate_func1,
        'column':'proteinID',
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'goldstandard/express')

    },
    {
        'id':9604717,
        'scoring_func':None,
        'validation_func':validate_func2_3,
        'column':'proteinID',
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'goldstandard/express/prospective_ova_pro_gold_express.txt')
    },
    {
        'id':9604718,
        'scoring_func':None,
        'validation_func':validate_func2_3,
        'column':'phosphoID',
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'goldstandard/express/prospective_ova_phospho_gold_express.txt')

    },
# Proteogenomics Subchallenge 1 Internal (9606530)
# Proteogenomics Subchallenge 2 Internal (9606531)
# Proteogenomics Subchallenge 3 Internal (9606532)
    {
        'id':9606530,
        'scoring_func':score1,
        'validation_func':validate_func1,
        'column':'proteinID',
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'goldstandard/internal')

    },
    {
        'id':9606531,
        'scoring_func':score2,
        'validation_func':validate_func2_3,
        'column':'proteinID',
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'goldstandard')
    },
    {
        'id':9606532,
        'scoring_func':score3,
        'validation_func':validate_func2_3,
        'column':'phosphoID',
        'goldstandard_path':os.path.join(os.path.dirname(os.path.abspath(__file__)),'goldstandard')
    }
]
evaluation_queue_by_id = {q['id']:q for q in evaluation_queues}


## define the default set of columns that will make up the leaderboard
LEADERBOARD_COLUMNS = [
    dict(name='objectId',      display_name='ID',      columnType='STRING', maximumSize=20),
    dict(name='userId',        display_name='User',    columnType='STRING', maximumSize=20, renderer='userid'),
    dict(name='entityId',      display_name='Entity',  columnType='STRING', maximumSize=20, renderer='synapseid'),
    dict(name='versionNumber', display_name='Version', columnType='INTEGER'),
    dict(name='name',          display_name='Name',    columnType='STRING', maximumSize=240),
    dict(name='team',          display_name='Team',    columnType='STRING', maximumSize=240)]

## Here we're adding columns for the output of our scoring functions, score,
## rmse and auc to the basic leaderboard information. In general, different
## questions would typically have different scoring metrics.
leaderboard_columns = {}
for q in evaluation_queues:
    leaderboard_columns[q['id']] = LEADERBOARD_COLUMNS + [
        dict(name='score',         display_name='Score',   columnType='DOUBLE'),
        dict(name='rmse',          display_name='RMSE',    columnType='DOUBLE'),
        dict(name='auc',           display_name='AUC',     columnType='DOUBLE')]

## map each evaluation queues to the synapse ID of a table object
## where the table holds a leaderboard for that question
leaderboard_tables = {}


def validate_submission(syn, evaluation, submission):
    """
    Find the right validation function and validate the submission.

    :returns: (True, message) if validated, (False, message) if
              validation fails or throws exception
    """
    config = evaluation_queue_by_id[int(evaluation.id)]
    validation_func = config['validation_func']

    status = syn.getSubmissionStatus(submission)
    #Filter through to fiind PREDICTION_FILE
    synId = filter(lambda input: input.get('key', None) == "PREDICTION_FILE", status.annotations['stringAnnos'])[0]
    entity = syn.get(synId['value'])

    dirname = os.path.dirname(entity.path)
    zfile = zipfile.ZipFile(entity.path)

    for name in zfile.namelist():
      zfile.extract(name, dirname)

    results, validation_message = validation_func(dirname, config['goldstandard_path'], config['column'])

    return True, validation_message


def score_submission(syn, evaluation, submission):
    """
    Find the right scoring function and score the submission

    :returns: (score, message) where score is a dict of stats and message
              is text for display to user
    """
    config = evaluation_queue_by_id[int(evaluation.id)]
    scoring_func = config['scoring_func']

    status = syn.getSubmissionStatus(submission)
    #Filter through to fiind PREDICTION_FILE
    synId = filter(lambda input: input.get('key', None) == "PREDICTION_FILE", status.annotations['stringAnnos'])[0]
    entity = syn.get(synId['value'])

    dirname = entity.cacheDir
    zfile = zipfile.ZipFile(entity.path)

    for name in zfile.namelist():
      zfile.extract(name, dirname)


    if scoring_func is not None:
        pred = scoring_func(dirname,config['goldstandard_path'])
    #Make sure to round results to 3 or 4 digits
        #return(dict(corr=round(corr,4), rmse=round(nrmse,4)), "You submission was scored.\ncorr: %s\nnrmse: %s" %(round(corr,4),round(nrmse,4)))
        return(pred, "Your submission was scored. You will receive your scores at the end of the round")
    else:
        return(dict(), "Your prediction file is in the correct format and can be scored.  Please feel free to submit to the real challenge queues.")



