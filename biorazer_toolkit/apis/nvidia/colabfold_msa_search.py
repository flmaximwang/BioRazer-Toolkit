#!/usr/bin/env python3
import asyncio
import json
from typing import Any, Dict, Optional
from fastapi import HTTPException
import httpx
import os
from pathlib import Path
from enum import StrEnum
import logging
from concurrent.futures import ThreadPoolExecutor

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


STATUS_URL = "https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/{task_id}"
PUBLIC_URL = "https://health.api.nvidia.com/v1/biology/colabfold/msa-search/predict"


def redirect_log(log: str, mode="a"):
    """
    Redirect log messages to a file.
    """
    file_handler = logging.FileHandler(log, mode=mode)
    # 清除现有的处理器，避免重复日志
    if logger.hasHandlers():
        logger.handlers.clear()
    logger.addHandler(file_handler)


def acquire_key() -> str:
    """
    Acquire the NVCF Run Key from the environment.
    """
    if os.environ.get("NVCF_RUN_KEY", None) is None:
        raise Exception("Error: Must set NVCF_RUN_KEY environment variable.")
    return os.environ.get("NVCF_RUN_KEY")


async def make_nvcf_call(
    function_url: str,
    data: Dict[str, Any],
    additional_headers: Optional[Dict[str, Any]] = None,
    NVCF_POLL_SECONDS: int = 10,
    MANUAL_TIMEOUT_SECONDS: int = 20,
) -> Dict:
    """
    Make a call to NVIDIA Cloud Functions using long-polling,
    which allows the request to patiently wait if there are many requests in the queue.
    """
    async with httpx.AsyncClient() as client:
        headers = {
            "Authorization": f"Bearer {acquire_key()}",
            "NVCF-POLL-SECONDS": f"{NVCF_POLL_SECONDS}",
            "Content-Type": "application/json",
        }
        if additional_headers is not None:
            headers.update(additional_headers)
        logger.debug(
            f"Headers: {dict(**{h: v for h, v  in headers.items() if 'Authorization' not in h})}"
        )
        # TIMEOUT must be greater than NVCF-POLL-SECONDS
        logger.debug(f"Making NVCF call to {function_url}")
        logger.debug(f"Data: {data}")
        response = await client.post(
            function_url, json=data, headers=headers, timeout=MANUAL_TIMEOUT_SECONDS
        )
        logger.debug(f"NVCF response: {response.status_code, response.headers}")

        if response.status_code == 202:
            # Handle 202 Accepted response
            task_id = response.headers.get("nvcf-reqid")
            while True:
                ## Should return in 5 seconds, but we set a manual timeout in 10 just in case
                status_response = await client.get(
                    STATUS_URL.format(task_id=task_id),
                    headers=headers,
                    timeout=MANUAL_TIMEOUT_SECONDS,
                )
                if status_response.status_code == 200:
                    return status_response.status_code, status_response
                elif status_response.status_code in [400, 401, 404, 422, 500]:
                    raise HTTPException(
                        status_response.status_code,
                        "Error while waiting for function:\n",
                        response.text,
                    )
        elif response.status_code == 200:
            return response.status_code, response
        else:
            raise HTTPException(status_code=response.status_code, detail=response.text)


async def call(
    sequence: str,
    e_value=0.0001,
    iterations=1,
    databases=["Uniref30_2302", "PDB70_220313", "colabfold_envdb_202108"],
    output_json="colabfold_msa_search.json",
    api_key: str = None,
):
    """
    async function to call the NVCF ColabFold MSA Search API. By default, results will be saved to output_json if provided.
    """

    if api_key is None:
        raise (Exception("Error: Must provide api_key argument."))
    os.environ["NVCF_RUN_KEY"] = api_key

    # Initial request
    ## Note: headers are set in make_nvcf_call function
    data = {
        "sequence": sequence,
        "e_value": e_value,
        "iterations": iterations,
        "databases": databases,
        "output_alignment_formats": ["a3m", "fasta"],
    }
    logger.info("Making request to NVCF...")
    code, response = await make_nvcf_call(function_url=PUBLIC_URL, data=data)

    if code == 200:
        logger.info(f"Request succeeded, returned {code}")
        response_dict = response.json()
        if not output_json is None:
            output_path = Path(output_json)
            with open(output_path, "w") as f:
                json.dump(response_dict, f, indent=4)
            logger.info(f"Wrote output to {output_path.resolve()}")
        ## print the dictionaries in the alignments portion of the response:
        logger.info(
            f"The returned databases were: {list(response_dict['alignments'].keys())} ."
        )
        for database in databases:
            if database not in response_dict["alignments"]:
                logger.warning(f"Database {database} not found in response.")
                continue
            else:
                ## print the file formats returned:
                logger.info(
                    f"The returned formats were: {list(response_dict['alignments']['Uniref30_2302'].keys())} ."
                )
                ## print the length of the FASTA-formatted alignment:
                logger.info(
                    f"The returned FASTA contained {len(response_dict['alignments']['Uniref30_2302']['fasta']['alignment'])} characters."
                )

        del os.environ["NVCF_RUN_KEY"]
        return response_dict
    else:
        logger.error(f"Request failed with status code {code}")
        del os.environ["NVCF_RUN_KEY"]
        raise HTTPException(status_code=code, detail=response.text)


def run(
    sequence: str,
    e_value=0.0001,
    iterations=1,
    databases=["Uniref30_2302", "PDB70_220313", "colabfold_envdb_202108"],
    api_key: str = None,
    output_json="colabfold_msa_search.json",
):
    """
    Synchronous wrapper to call the NVCF ColabFold MSA Search API. Use this function in .py scripts. Use the async call in Jupyter notebooks.
    """

    res = asyncio.run(
        call(
            sequence=sequence,
            e_value=e_value,
            iterations=iterations,
            databases=databases,
            output_json=output_json,
            api_key=api_key,
        )
    )

    return res
