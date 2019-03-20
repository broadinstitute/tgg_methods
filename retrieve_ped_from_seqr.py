import getpass
import json
import requests
from pprint import pprint

SEQR_URL = "https://seqr.broadinstitute.org"
LOGIN_URL = SEQR_URL+"/login"
INDIVIDUAL_API_URL = SEQR_URL+"/api/project/<project_guid>/export_project_individuals?file_format=json"

seqr_username = "mwilson@broadinstitute.org"
seqr_password = getpass.getpass()     # prompt for password

s = requests.Session()
response = s.get(LOGIN_URL)    # get CSRF cookies
response.raise_for_status()

csrf_token = s.cookies.get('csrftoken')
response = s.post(
  LOGIN_URL,
  data={
    'username_or_email': seqr_username,
    'password': seqr_password,
    'csrfmiddlewaretoken': csrf_token
  },
  headers={
    'Content-Type': 'application/x-www-form-urlencoded'
  })    # submit login

response.raise_for_status()

project_guid = input("Enter seqr project guid (eg. R0390_1000_genomes_demo): ")

response = s.get(INDIVIDUAL_API_URL.format(**locals()))   # get json from API
response.raise_for_status()
pprint(response.text)
# parse json response
#json_data = [json.loads(row) for row in response.content.strip().split("\n")]
#pprint(json_data)
