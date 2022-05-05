import argparse
import logging

import requests
from typing import Set

logging.basicConfig(
    format="%(asctime)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

SEQR_URL = "https://seqr.broadinstitute.org"
GET_FAMILIES_URL = "get_families"
GET_INDIVIDUALS_URL = "get_individuals"


def pull_project_peds(session_id: str, projects: Set[str]):
    """
    Use the seqr API to retrieve a list of projects' pedigrees. 

    Writes a single pedigree for all projects ('seqr_pedigrees.txt'), using the retrieved pedigrees. 
    Failing projects are written to 'projects_not_pulled.txt'.

    :param session_id: sessionid cookie used for seqr login
    :param projects:  Set of seqr project GUIDs
    """
    with requests.Session() as s:
        # Initialize session with authenticated cookie
        s.cookies.set_cookie(requests.cookies.create_cookie("sessionid", session_id))
        with open("seqr_pedigrees.txt", "w") as final_ped, open(
            "projects_not_pulled.txt", "w"
        ) as errors:
            final_ped.write(
                "Project_GUID\tFamily_ID\tIndividual_ID\tPaternal_ID\tMaternal_ID\tSex\n"
            )
            for project_guid in projects:
                project_guid = project_guid.rstrip()
                family_r = s.get(
                    f"{SEQR_URL}/api/project/{project_guid}/{GET_FAMILIES_URL}"
                )
                individual_r = s.get(
                    f"{SEQR_URL}/api/project/{project_guid}/{GET_INDIVIDUALS_URL}"
                )
                if (family_r.status_code != requests.codes["ok"]) or (
                    individual_r.status_code != requests.codes["ok"]
                ):
                    if family_r.status_code != requests.codes["ok"]:
                        logger.info(
                            f"Could not get family IDs for {project_guid} in seqr."
                        )
                    if individual_r.status_code != requests.codes["ok"]:
                        logger.info(
                            f"Could not get individual IDs for {project_guid} in seqr."
                        )
                    errors.write(f"{project_guid}\n")
                else:
                    # `family_r` contains a high level field `familiesByGuid`,
                    # which contains fields for individual families.
                    # The fields within `familiesByGuid` also contain various fields (e.g., analysisStatus, analyzedBy, etc).
                    # Importantly, familiesByGuid has the fields `familyGuid` (keyword type family name) 
                    # and `familyId`, which is the family name used in seqr
                    fams = family_r.json()["familiesByGuid"]
                    fams_dict = {}
                    for fam in fams.values():
                        fams_dict[fam["familyGuid"]] = fam["familyId"]
                    # `individual_r` contains the high level field `individualsByGuid`, which contains fields for individuals
                    # such as `paternalId`, `maternalId`, and `familyGuid`. 
                    # Critically, it does not contain `familyId`, 
                    # which is why we need to build a mapping from `familyGuid` to `familyId`.
                    inds = individual_r.json()["individualsByGuid"]
                    for ind in inds.values():
                        family_id = fams_dict[ind["familyGuid"]]
                        final_ped.write(
                            "\t".join(
                                map(
                                    str,
                                    [
                                        ind["projectGuid"],
                                        family_id,
                                        ind["individualId"],
                                        ind["paternalId"],
                                        ind["maternalId"],
                                        ind["sex"] if ind["sex"] != "U" else "None",
                                    ],
                                )
                            )
                            + "\n"
                        )


def main(args):
    session_id = args.session_id
    projects = set([])

    with open(args.projects, "r") as project_list:
        logger.info("Retreiving the following projects pedigrees:")
        for project_guid in project_list:
            project_guid = project_guid.rstrip()
            logger.info(f"{project_guid}")
            projects.add(project_guid)

    pull_project_peds(session_id, projects)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--session-id",
        help="Session ID cookie retrieved from a logged in seqr session",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--projects",
        help="Text file with a list of seqr projects to pull",
        required=True,
    )
    args = parser.parse_args()
    main(args)
