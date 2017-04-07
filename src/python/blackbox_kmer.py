#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

from collections import Counter

from sklearn import metrics
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB

from sklearn.model_selection import train_test_split
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler

from imblearn.ensemble import EasyEnsemble

from sklearn.externals import joblib

import os

class naive:
  def __init__(self):
      pass

  def fit(self, X, Y):
      pass
      
  def predict(self, X):
      xSize = X.size / 6
      Y = np.zeros((xSize, 1))
      for i in range(xSize):
        entry = X[i]
        countObservedBiasCorrected = entry[4]
        countExpectedPusm = entry[5]
        if (countObservedBiasCorrected < 0.5 * countExpectedPusm):
            Y[i][0] = 2 #UNTRUSTED
        elif (countObservedBiasCorrected > 1.5 * countExpectedPusm):
            Y[i][0] = 0 #REPEAT
        else:
            Y[i][0] = 1 #TRUSTED
      return Y
      
      
class statistical:
  def __init__(self):
      pass

  def fit(self, X, Y):
      pass
      
  def predict(self, X):
      xSize = X.size / 6
      Y = np.zeros((xSize, 1))
      for i in range(xSize):
        entry = X[i]
        zscore = entry[0]
        if (zscore < -2):
            Y[i][0] = 2 #UNTRUSTED
        elif (zscore > 2):
            Y[i][0] = 0 #REPEAT
        else:
            Y[i][0] = 1 #TRUSTED
        pass   
      return Y


class classifier:
  def __init__(self):
    self.features = []
    self.classes = []
    #self.models = [GaussianNB(), DecisionTreeClassifier(), DecisionTreeClassifier(class_weight = 'balanced'), RandomForestClassifier(), RandomForestClassifier(class_weight = 'balanced'), LogisticRegression(), LogisticRegression(class_weight = 'balanced')]#, AdaBoostClassifier(), AdaBoostClassifier(DecisionTreeClassifier(class_weight = 'balanced'))]
    #self.modelnames = ['GaussianNB', 'DecisionTreeClassifier', 'DecisionTreeClassifier(balanced)', 'RandomForestClassifier', 'RandomForestClassifier(balanced)', 'LogisticRegression', 'LogisticRegression(balanced)']#, 'AdaBoostClassifier', 'AdaBoostClassifier(balanced)']
    
    
    self.models = [naive(), statistical(), GaussianNB(), DecisionTreeClassifier(class_weight = 'balanced'), RandomForestClassifier(class_weight = 'balanced'), LogisticRegression(class_weight = 'balanced')]#, AdaBoostClassifier(), AdaBoostClassifier(DecisionTreeClassifier(class_weight = 'balanced'))]
    self.modelnames = ['naive', 'statistical', 'GaussianNB', 'DecisionTreeClassifier', 'RandomForestClassifier', 'LogisticRegression(balanced)']#, 'AdaBoostClassifier', 'AdaBoostClassifier(balanced)']
    
    
    self.best_model = naive()
  
  def add(self, num1, num2):
    print "miau"
    return num1 + num2

  def read_data(self, input_file):
    print("Reading data...")
    df = pd.read_csv(input_file, header = 0, delimiter= ';')
    X = df[self.features].values
    Y = df['type'].values
    return (X, Y)
    
  #takes std::vector<std::string>
  def set_features(self, features):
    self.features = features
    
  #takes std::vector<int>
  def set_classes(self, classes):
    self.classes = classes
    
  def balance_dataset(self, X, Y):
    X_new = X
    Y_new = Y
    overSampler = RandomOverSampler()
    underSampler = RandomUnderSampler()
    #sm = EasyEnsemble()
    #X_refit, Y_refit = sm.fit_sample(X, Y)
    #print('Resampled dataset shape {}'.format(Counter(Y_refit[0])))
    #X, Y = X_refit[0], Y_refit[0]
    classCounts = Counter(Y)
    print('Original training dataset shape {}'.format(classCounts))
    avg = 0
    minCount = classCounts[self.classes[0]]
    maxCount = classCounts[self.classes[0]]
    for i in self.classes:
        avg = avg + classCounts[i]
        if classCounts[i] < minCount:
            minCount = classCounts[i]
        if classCounts[i] > maxCount:
            maxCount = classCounts[i]
    avg = avg // len(classCounts)
    print("Rounded-down average class count in training dataset: " + str(avg))
    print("minCount: " + str(minCount))
    print("maxCount: " + str(maxCount))
    rate = avg / float(maxCount)
    print("rate: " + str(rate))
    underSampler = RandomUnderSampler(ratio = rate)
    X_new, Y_new = underSampler.fit_sample(X_new, Y_new)
    classCounts = Counter(Y_new)
    print('Class counts after undersampling {}'.format(classCounts))
    #avg = 0
    #minCount = classCounts[0]
    #maxCount = classCounts[0]
    #for i in range(len(classCounts)):
    #    avg = avg + classCounts[i]
    #    if classCounts[i] < minCount:
    #        minCount = classCounts[i]
    #    if classCounts[i] > maxCount:
    #        maxCount = classCounts[i]
    #avg = avg // len(classCounts)
    #rate = minCount / float(avg)
    #print("rate: " + str(rate))
    #overSampler = RandomOverSampler(ratio = rate)
    #print("I am here1")
    X_new, Y_new = overSampler.fit_sample(X_new, Y_new)
    #print("I am here2")
    return X_new, Y_new

  def choose_best_classifier(self, X_train, X_test, y_train, y_test, bestF, useSampling):
    if useSampling == True:
        print("Trying to choose the best classifier with presampling the training dataset...")
        X_train, y_train = self.balance_dataset(X_train, y_train)
    else:
        print("Trying to choose the best classifier without presampling the training dataset...")
    expected = y_test
    best_fscore = bestF
    best_num = 0
    
    for i in range(len(self.models)):
        # Fit the model
        print("Fitting model: " + self.modelnames[i] + "...")
        self.models[i].fit(X_train, y_train)
        # Evaluate the model
        print("Predicting class labels with model: " + self.modelnames[i] + "...")
        predicted = self.models[i].predict(X_test)
        # print some metrics
        print(metrics.confusion_matrix(expected, predicted))
        print(metrics.classification_report(expected, predicted))
        fscore = metrics.f1_score(expected, predicted, average = 'macro')
        print("The average F-Score for model " + self.modelnames[i] + " with equal class weights is: " + str(fscore))
        if fscore > best_fscore:
            self.best_model = self.models[i]
            best_fscore = fscore
            best_num = i
    print("And the winner is: " + self.modelnames[best_num])
    return best_fscore
      
  #takes std::string
  # TODO: Mix oversampling and undersampling in order to get (original average count) many entries per class
  def set_csv_file(self, csv_path):
    print(csv_path)
    (X, Y) = self.read_data(csv_path)
    print('Original dataset shape {}'.format(Counter(Y)))
    # Split into training and testing data
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.34)
    bestF = 0
    # Do not prebalance the dataset before training, use sample weights
    bestF = self.choose_best_classifier(X_train, X_test, y_train, y_test, bestF, False)
    #bestFOld = bestF
    # Try to prebalance the dataset before training
    #bestF = self.choose_best_classifier(X_train, X_test, y_train, y_test, bestF, True)
    #if bestF > bestFOld:
    #    print("Resampling was good.")
    #else:
    #    print("Resampling was useless or even harmful.")
    
  #takes std::string
  def store_classifier(self, filename):
      joblib.dump(self.best_model, filename, compress=9)
      print("Stored classifier.")
      
  #takes std::string
  def load_classifier(self, filename):
      self.best_model = joblib.load(filename)
      print("Loaded classifier.")

  #takes std::vector<double> of size len(_features), returns std::vector<double> of size len(_classes)
  def proba(self, feature_vector):
    matrix = []
    matrix.append(feature_vector)
    probs = self.best_model.predict_proba(matrix)[0]
    res = [float(i) for i in probs]
    return res

  #takes std::vector<double>, returns a single int
  def classify(self, feature_vector):
    matrix = []
    matrix.append(feature_vector)
    return self.best_model.predict(matrix)[0]
