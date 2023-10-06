from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import GroupKFold

from sklearn.ensemble import RandomForestClassifier

def LR_withGridSearchCV(X_train,y_train,n_jobs=20):
    LRtopDCAs_model = LogisticRegression(penalty="l1",solver='liblinear', random_state=0,
                                         max_iter=100,
                                 ) # parameter from paper Large-scale discovery of protein interactions at residue resolution using co-evolution calculated from genomic sequences
    # alreay tried parameters 
    # 'C':[0.1,1, 10]
    # 'class_weight':[{1: w} for w in [0.1, 1, 10,50,100]]
    parameters = {'C':[0.1,1, 10],
                 'class_weight':[{1: w} for w in [0.1, 1, 10,50,100]]}
    clrtopDCAs = GridSearchCV(LRtopDCAs_model, parameters, cv=5,
                         n_jobs=n_jobs,
                      #scoring=sklearn.metrics.make_scorer(sklearn.metrics.f1_score)
                             )
    clrtopDCAs.fit(X_train,y_train)
    return(clrtopDCAs)

def RF_withGridSearchCV(X_train,y_train,n_jobs=20):
    RFtopDCAs_model = RandomForestClassifier(random_state=0)
    random_grid = {'n_estimators': [10,100,500],
                   #'max_features': max_features,
                   #'max_depth': max_depth,
                   'min_samples_split': [2,5,10],
                   'min_samples_leaf': [1,2,4],
                   #'bootstrap': bootstrap
                   # 'class_weight':[{1: w} for w in [None,"balanced","balanced_subsample"]]
                   # 'class_weight':[None,"balanced","balanced_subsample"],
                  }
    
    cRFtopDCAs = GridSearchCV(RFtopDCAs_model, random_grid, cv=5,
                         n_jobs=n_jobs,
                      #scoring=sklearn.metrics.make_scorer(sklearn.metrics.f1_score)
                             )
    cRFtopDCAs.fit(X_train,y_train)
    return(cRFtopDCAs)


def LR_withGridSearchCV_GroupKFold(X_train,y_train,group,n_jobs=20,
                                  scoring_metrics=None):    
    import warnings
    from sklearn.exceptions import ConvergenceWarning
    #https://stackoverflow.com/questions/15933741/how-do-i-catch-a-numpy-warning-like-its-an-exception-not-just-for-testing/15934081#15934081
    # https://stackoverflow.com/questions/48100939/how-to-detect-a-scikit-learn-warning-programmatically
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            LRtopDCAs_model = LogisticRegression(penalty="l1",solver='liblinear', random_state=0,
                                                 max_iter=100,
                                                
                                         ) # parameter from paper Large-scale discovery of protein interactions at residue resolution using co-evolution calculated from genomic sequences
            # alreay tried parameters 
            # 'C':[0.1,1, 10]
            # 'class_weight':[{1: w} for w in [0.1, 1, 10,50,100]]
            parameters = {'C':[0.1,1, 10],
                         'class_weight':[{1: w} for w in [0.1, 1, 10,50,100]]}
            gkf = GroupKFold(n_splits=5).split(X_train,y_train,group)
            clrtopDCAs = GridSearchCV(LRtopDCAs_model, parameters, cv=gkf,
                                 n_jobs=n_jobs,
                                scoring=scoring_metrics
                              #scoring=sklearn.metrics.make_scorer(sklearn.metrics.f1_score)
                                     )
            clrtopDCAs.fit(X_train,y_train)
        except ConvergenceWarning as w: 
            print("LR ConvergenceWarning: ",w)
            LRtopDCAs_model = LogisticRegression(penalty="l1",solver='liblinear', random_state=0,
                                                 max_iter=500,
                                         ) 
            parameters = {'C':[0.1,1, 10],
                         'class_weight':[{1: w} for w in [0.1, 1, 10,50,100]]}
            gkf = GroupKFold(n_splits=5).split(X_train,y_train,group)
            clrtopDCAs = GridSearchCV(LRtopDCAs_model, parameters, cv=gkf,
                                 n_jobs=n_jobs,
                                scoring=scoring_metrics
                              #scoring=sklearn.metrics.make_scorer(sklearn.metrics.f1_score)
                                     )
            clrtopDCAs.fit(X_train,y_train)
            # assert 1==2

    return(clrtopDCAs)

def RF_withGridSearchCV_GroupKFold(X_train,y_train,group,n_jobs=20,
                                  scoring_metrics=None):
    RFtopDCAs_model = RandomForestClassifier(random_state=0,n_jobs=4)
    random_grid = {'n_estimators': [10,100,500],
                   'max_features': [2,3,5],#max_features,
                   'max_depth':[5,8,8], #  ,10,max_depth, 9,10,12,15, seem alredy lead overfitting 
                   # 'min_samples_split': [2,5,10],
                   # 'min_samples_leaf': [1,2,4],
                    'class_weight':[{1: w} for w in [1, 10,50,]],
                   #'class_weight':[{0: w} for w in [1, 10,50,100]],
# 'class_weight':[{0: w} for w in [None,"balanced","balanced_subsample", 
#                                 10,50,100]],
                   #'bootstrap': bootstrap
                  }
    gkf = GroupKFold(n_splits=5).split(X_train,y_train,group)
    cRFtopDCAs = GridSearchCV(RFtopDCAs_model, random_grid, cv=gkf,
                         n_jobs=n_jobs,
                        scoring=scoring_metrics
                      #scoring=sklearn.metrics.make_scorer(sklearn.metrics.f1_score)
                             )
    cRFtopDCAs.fit(X_train,y_train)
    return(cRFtopDCAs)

