"""Models compared in this repo. Hyperparameters are fixed, not tuned on test folds."""

from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def logistic():
    """Logistic regression: the one-feature pLDDT baselines, and a linear check on all features."""
    return make_pipeline(
        SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(max_iter=1000)
    )


def gradient_boosting(seed: int = 0):
    """Shallow, regularised boosting (small data), Platt-calibrated on inner folds."""
    gb = HistGradientBoostingClassifier(
        max_depth=3,
        learning_rate=0.05,
        max_iter=200,
        min_samples_leaf=20,
        l2_regularization=1.0,
        random_state=seed,
    )
    return CalibratedClassifierCV(gb, method="sigmoid", cv=3)
