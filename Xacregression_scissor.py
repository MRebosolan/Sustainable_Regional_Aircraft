import pandas as pd
from sklearn import linear_model
import statsmodels.api as sm

def Xacregression(BA, taper, sweepbeta):
    Stock_Market = {'BA': [6, 6, 6, 4, 4, 4],
                    'taper': [0.2, 0.3, 0.2, 0.2, 0.3, 0.2],
                    'sweepbeta': [50, 50, 40, 50, 50, 40],        
                    'xac/c':[0.45, 0.41, 0.38, 0.4, 0.37, 0.35]}
    
    df = pd.DataFrame(Stock_Market,columns=['BA', 'taper', 'sweepbeta', 'xac/c'])
    
    X = df[['BA', 'taper', 'sweepbeta']] # here we have 2 variables for multiple regression. If you just want to use one variable for simple linear regression, then use X = df['Interest_Rate'] for example.Alternatively, you may add additional variables within the brackets
    Y = df['xac/c']
     
    # with sklearn
    regr = linear_model.LinearRegression()
    regr.fit(X, Y)
    
    # print('Intercept: \n', regr.intercept_)
    # print('Coefficients: \n', regr.coef_)
    
    # prediction with sklearn
    # BA = 5.1
    # taper = 0.2
    # sweepbeta = 40
    # print ('Predicted Xac/c: \n', regr.predict([[BA, taper, sweepbeta]]))
    if BA>7 or BA<3:
        return print('BA is out of range, keep it between 7 and 3')
    elif taper>0.4 or taper <0.125:
        return print('taper is out of range for xac estimation')
    elif sweepbeta >60 or sweepbeta < 30:
        return print('sweepbeta is out of range for xac estimation')
    else:
        return regr.predict([[BA, taper, sweepbeta]])[0]

def Xacregression_app(BA, taper, sweepbeta):
    Stock_Market = {'BA': [6, 6, 6, 8, 8, 8],
                    'taper': [0.2, 0.3, 0.2, 0.2, 0.3, 0.2],
                    'sweepbeta': [30, 30, 20, 30, 30, 20],        
                    'xac/c':[0.325, 0.305, 0.28, 0.35, 0.325, 0.3]}
    
    df = pd.DataFrame(Stock_Market,columns=['BA', 'taper', 'sweepbeta', 'xac/c'])
    
    X = df[['BA', 'taper', 'sweepbeta']] # here we have 2 variables for multiple regression. If you just want to use one variable for simple linear regression, then use X = df['Interest_Rate'] for example.Alternatively, you may add additional variables within the brackets
    Y = df['xac/c']
     
    # with sklearn
    regr = linear_model.LinearRegression()
    regr.fit(X, Y)
    
    # print('Intercept: \n', regr.intercept_)
    # print('Coefficients: \n', regr.coef_)
    
    # prediction with sklearn
    # BA = 5.1
    # taper = 0.2
    # sweepbeta = 40
    # print ('Predicted Xac/c: \n', regr.predict([[BA, taper, sweepbeta]]))
    if BA>9 or BA<5:
        return print('BA is out of range, keep it between 7 and 3')
    elif taper>0.4 or taper <0.125:
        return print('taper is out of range for xac estimation')
    elif sweepbeta >40 or sweepbeta < 15:
        return print('sweepbeta is out of range for xac estimation')
    else:
        return regr.predict([[BA, taper, sweepbeta]])[0]
    
    
def moment_flap_extension(BA, taper, sweepbeta):
    Stock_Market = {'BA': [6, 6, 6, 8, 8, 8],
                    'taper': [0.2, 0.3, 0.2, 0.2, 0.3, 0.2],
                    'sweepbeta': [30, 30, 20, 30, 30, 20],        
                    'xac/c':[0.325, 0.305, 0.28, 0.35, 0.325, 0.3]}
    
    df = pd.DataFrame(Stock_Market,columns=['BA', 'taper', 'sweepbeta', 'xac/c'])
    
    X = df[['BA', 'taper', 'sweepbeta']] # here we have 2 variables for multiple regression. If you just want to use one variable for simple linear regression, then use X = df['Interest_Rate'] for example.Alternatively, you may add additional variables within the brackets
    Y = df['xac/c']
     
    # with sklearn
    regr = linear_model.LinearRegression()
    regr.fit(X, Y)
    
    # print('Intercept: \n', regr.intercept_)
    # print('Coefficients: \n', regr.coef_)
    
    # prediction with sklearn
    # BA = 5.1
    # taper = 0.2
    # sweepbeta = 40
    # print ('Predicted Xac/c: \n', regr.predict([[BA, taper, sweepbeta]]))
    if BA>9 or BA<5:
        return print('BA is out of range, keep it between 7 and 3')
    elif taper>0.4 or taper <0.125:
        return print('taper is out of range for xac estimation')
    elif sweepbeta >40 or sweepbeta < 15:
        return print('sweepbeta is out of range for xac estimation')
    else:
        return regr.predict([[BA, taper, sweepbeta]])[0]