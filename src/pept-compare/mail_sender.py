import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import os

me  = 'peptimetric@gmail.com'
recipient = 'peptimetric@gmail.com'
subject = 'Feedback from peptimetric webapp'
password = 'P_metric9998'

email_server_host = 'smtp.gmail.com'
port = 587
email_username = me
email_password = password

msg = MIMEMultipart('alternative')
msg['From'] = me
msg['To'] = recipient
msg['Subject'] = subject

def send_email(email_body):
    msg.attach(MIMEText(email_body, 'html'))
    server = smtplib.SMTP(email_server_host, port)
    server.ehlo()
    server.starttls()
    server.login(email_username, email_password)
    server.sendmail(me, recipient, msg.as_string())
    server.close()