from flask_wtf.form import FlaskForm
from wtforms.fields import EmailField
from wtforms.fields.simple import TextAreaField
from wtforms import validators

class SupportForm(FlaskForm):
    email = EmailField('Email',
                       [validators.DataRequired(), validators.Email()])
    body = TextAreaField('Body')