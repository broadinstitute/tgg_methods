import requests
import getpass
import argparse


def pull_project_peds(args):
    """
    Uses the seqr API to retrieve a list of projects' pedigrees. Writes a single pedigree,
    'pulled_seqr_pedigrees.txt', using the retrieved pedigrees. Failing projects are written
    to 'projects_not_pulled.txt'
    """
    email=args.email
    projects=args.projects

    with requests.Session() as s:
        p = s.post(
            "https://seqr.broadinstitute.org/api/login", 
            json={
                'email': email,
                'password': getpass.getpass(),
            }
        )
        final_ped = open('seqr_pedigrees.txt', 'w')
        final_ped.write('Project_GUID\tFamily_ID\tIndividual_ID\tPaternal_ID\tMaternal_ID\tSex\n') 
        errors = open('projects_not_pulled.txt', 'w')

        with open(projects, 'r') as project_list:
            for project_guid in project_list:
                project_guid = project_guid.rstrip()
                r = s.get(f'https://seqr.broadinstitute.org/api/project/{project_guid}/details') 
                if r.status_code != requests.codes['ok']:
                    errors.write(f'{project_guid}\n')
                else:
                    fams = r.json()['familiesByGuid']
                    fams_dict = {}
                    for key in fams:
                        fams_dict[fams[key]['familyGuid']] = fams[key]['familyId']
                    inds = r.json()['individualsByGuid']
                    for key in inds:
                        family_id = fams_dict[inds[key]['familyGuid']]
                        project_guid = inds[key]['projectGuid']
                        individual_id = inds[key]['individualId']
                        paternal_id = inds[key]['paternalId']
                        maternal_id = inds[key]['maternalId']
                        sex = inds[key]['sex']
                        final_ped.write(f'{project_guid}\t{family_id}\t{individual_id}\t{paternal_id}\t{maternal_id}\t{sex}\n')
        final_ped.close()
        errors.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--email', help="Email used for seqr login", required=True)
    parser.add_argument('-p', '--projects', help="Text file with a list of projects to pull", required=True)
    args = parser.parse_args()
    pull_project_peds(args)
    