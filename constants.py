from typing import List

class Predictor:
    def __init__(self, short_name, full_name):
        self.short_name = short_name
        self.full_name = full_name

    def __str__(self):
        return self.short_name
    
    def to_dict(self):
        return {"short_name": self.short_name, "full_name": self.full_name}

    def __repr__(self):
        return f"Predictor(short_name='{self.short_name}', full_name='{self.full_name}')"

    @classmethod
    def from_dict(cls, data):
        return cls(data['short_name'], data['full_name'])
    
class Class_One_Predictors:
    MixMHCpred = Predictor("MixMHCpred", "mixMHCpred 3.0")
    MHCflurry = Predictor("MHCflurry", "MHCflurry 2.0")
    NetMHCpan = Predictor("NetMHCpan", "netMHCpan 4.1 b")


class Class_Two_Predictors:
    MixMHC2pred = Predictor("MixMHC2pred", "MixMHC2pred-2.0")
    NetMHCpanII = Predictor("NetMHCpanII", "netMHCIIpan 4.3 e")

def get_all_predictors() -> List[Predictor]:
    """Fetch all predictors dynamically from Class_One_Predictors and Class_Two_Predictors."""
    from constants import Class_One_Predictors, Class_Two_Predictors  # Avoid circular imports
    return [
        getattr(Class_One_Predictors, attr)
        for attr in dir(Class_One_Predictors)
        if isinstance(getattr(Class_One_Predictors, attr), Predictor)
    ] + [
        getattr(Class_Two_Predictors, attr)
        for attr in dir(Class_Two_Predictors)
        if isinstance(getattr(Class_Two_Predictors, attr), Predictor)
    ]

class MHC_Class:
    One = "I"
    Two = "II"