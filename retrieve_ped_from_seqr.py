#!/usr/bin/env python
import getpass
import requests
import argparse

# Hack for 2.X and 3.X to use input()
try:
    input = raw_input
except NameError:
    pass


def pull_project_peds(args):
    '''
    Uses the seqr API to retrieve a list of projects' pedigrees. Writes a single pedigree,
    'compiled_pedigree.txt', using the retrieved pedigrees.
    :param args: List of projects to be retrieved
    '''
    projects = args.projects
    login_url = "https://seqr.broadinstitute.org/login"
    individual_api_url = "https://seqr.broadinstitute.org/api/project/{project_guid}/export_project_individuals?file_format=tsv"

    seqr_username = 'mwilson@broadinstitute.org'#input("What is your seqr username? ")
    seqr_password = getpass.getpass()

    s = requests.Session()
    response = s.get(login_url, verify=False)  # get CSRF cookies
    response.raise_for_status()

    csrf_token = s.cookies.get('csrftoken')
    response = s.post(
        login_url,
        data={
            'username_or_email': seqr_username,
            'password': seqr_password,
            'csrfmiddlewaretoken': csrf_token
        },
        headers={
            'Content-Type': 'application/x-www-form-urlencoded'
        })  # submit login
    response.raise_for_status()

    temp_ped = open('temp_pedigree.txt', 'w')
    errors = open('projects_not_retrieved.txt', 'w')

    with open(projects, 'r') as project_list:
        for project_guid in project_list:
            project_guid = project_guid.rstrip()
            response = s.get(individual_api_url.format(**locals()))  # get tsv from API
            if response.status_code == requests.codes.ok:
                temp_ped.write(response.text)
            else:
                errors.write('{project_guid}\n'.format(**locals()))

    errors.close()
    temp_ped.close()

    final_ped = open('compiled_pedigree.ped', 'w')

    with open('temp_pedigree.txt', 'r') as ped:
        header = ped.readline()
        final_ped.write(header)
        for line in ped:
            if line != header:
                final_ped.write(line)

    final_ped.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--projects', help="File of projects in to retrieve, separated by newline", required=True)
    args = parser.parse_args()
    pull_project_peds(args)
