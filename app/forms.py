from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, SubmitField
from wtforms.validators import DataRequired, Email, Length, EqualTo
from flask_wtf.file import FileField, FileRequired

class InitialiserForm(FlaskForm):
    col_name = StringField("Peptide Column", validators=[DataRequired(),Length(min=2,max=55)])
    file_input = FileField("Peptide File", validators=[FileRequired()])
    submit = SubmitField("Register Now")