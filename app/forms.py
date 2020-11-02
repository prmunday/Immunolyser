from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, SubmitField, MultipleFileField
from wtforms.validators import DataRequired, Email, Length, EqualTo
from flask_wtf.file import FileField, FileRequired, DataRequired

class InitialiserForm(FlaskForm):
    sample_one_name = StringField("Sample One", validators=[DataRequired(),Length(min=2,max=55)])
    sample_one = MultipleFileField("Peptide File 1", validators=[DataRequired()])
    sample_two_name = StringField("Sample Two", validators=[DataRequired(),Length(min=2,max=55)])
    sample_two = MultipleFileField("Peptide File 2" , validators=[DataRequired()])
    submit = SubmitField("Submit")