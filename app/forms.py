from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, SubmitField, MultipleFileField, FieldList, FormField, IntegerField
from wtforms.validators import DataRequired, Email, Length, EqualTo
from flask_wtf.file import FileField, FileRequired, DataRequired

class InitialiserForm(FlaskForm):
    sample_name = StringField("Sample", validators=[DataRequired(),Length(min=2,max=55)])
    sample = MultipleFileField("Peptide File", validators=[DataRequired()])

class ParentForm(FlaskForm):

    children = FieldList(FormField(InitialiserForm), label='Children')
    add_sample = SubmitField(label='Add Sample')

    submit = SubmitField("Submit")
