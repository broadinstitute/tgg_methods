import argparse
import getpass
import logging

import requests

logging.basicConfig(
    format="%(asctime)s): %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def pull_project_peds(email: str, projects: set):
    """
    Uses the seqr API to retrieve a list of projects' pedigrees. Writes a single pedigree,
    'seqr_pedigrees.txt', using the retrieved pedigrees. Failing projects are written
    to 'projects_not_pulled.txt'
    :param email: email used for seqr login
    : param projects:  set of seqr project GUIDs
    """
    with requests.Session() as s:
        p = s.post(
            "https://seqr.broadinstitute.org/api/login",
            json={"email": email, "password": getpass.getpass(),},
        )
        with open("seqr_pedigrees.txt", "w") as final_ped, open("projects_not_pulled.txt", "w") as errors:
            final_ped.write(
                "Project_GUID\tFamily_ID\tIndividual_ID\tPaternal_ID\tMaternal_ID\tSex\n"
            )
            for project_guid in projects:
                project_guid = project_guid.rstrip()
                r = s.get(
                    f"https://seqr.broadinstitute.org/api/project/{project_guid}/details"
                )

                if r.status_code != requests.codes["ok"]:
                    logger.info(f"Could not find {project_guid} in seqr")
                    errors.write(f"{project_guid}\n")
                else:
                    fams = r.json()["familiesByGuid"]
                    fams_dict = {}
                    for fam in fams.values():
                        fams_dict[fam["familyGuid"]] = fam["familyId"]
                    inds = r.json()["individualsByGuid"]
                    for ind in inds.values():
                        family_id = fams_dict[ind["familyGuid"]]
                        project_guid = ind["projectGuid"]
                        individual_id = ind["individualId"]
                        paternal_id = ind["paternalId"]
                        maternal_id = ind["maternalId"]
                        sex = ind["sex"]
                        final_ped.write(
                            "\t".join(
                                [
                                    project_guid,
                                    family_id,
                                    individual_id,
                                    str(paternal_id),
                                    str(maternal_id),
                                    sex,
                                ]
                            )
                            + "\n"
                        )


def main(args):
    email = args.email
    projects = set([])

    with open(args.projects, "r") as project_list:
        logger.info("Retreiving the following projects pedigrees:")
        for project_guid in project_list:
            project_guid = project_guid.rstrip()
            logger.info(f"{project_guid}")
            projects.add(project_guid)

    pull_project_peds(email, projects)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--email", help="Email used for seqr login", required=True)
    parser.add_argument(
        "-p",
        "--projects",
        help="Text file with a list of seqr projects to pull",
        required=True,
    )
    args = parser.parse_args()
    main(args)
