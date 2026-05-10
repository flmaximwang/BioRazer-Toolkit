from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
import os
import tempfile
from typing import Callable, Iterator, TypeVar

from biotite.structure import AtomArray
from biorazer.structure.io.protein import STRUCT2CIF, STRUCT2PDB

T = TypeVar("T")

_WRITER_MAP = {
    "pdb": STRUCT2PDB,
    "cif": STRUCT2CIF,
}


def _normalize_file_format(file_format: str) -> str:
    normalized = file_format.lower().lstrip(".")
    if normalized not in _WRITER_MAP:
        supported = ", ".join(sorted(_WRITER_MAP))
        raise ValueError(
            f"Unsupported structure file format '{file_format}'. Supported formats: {supported}."
        )
    return normalized


def write_structure_file(
    atom_array: AtomArray,
    output_file: str | Path,
    file_format: str | None = None,
    **write_kwargs,
) -> Path:
    """
    Export an AtomArray into a structure file using biorazer exporters.

    Parameters
    ----------
    atom_array:
        Structure to export.
    output_file:
        Target output path.
    file_format:
        Explicit format, currently supports "pdb" and "cif".
        If not provided, infer from output file suffix.
    write_kwargs:
        Extra keyword arguments forwarded to biotite's set_structure.
    """
    output_path = Path(output_file)
    inferred_format = output_path.suffix or ".pdb"
    resolved_format = _normalize_file_format(file_format or inferred_format)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    writer_cls = _WRITER_MAP[resolved_format]
    writer_cls(output_file=output_path).write(atom_array, **write_kwargs)
    return output_path


@contextmanager
def atom_array_as_temp_file(
    atom_array: AtomArray,
    file_format: str = "pdb",
    temp_file_format: str | None = None,
    prefix: str = "biorazer_structure_",
    suffix: str | None = None,
    dir: str | Path | None = None,
    delete_file: bool = True,
    **write_kwargs,
) -> Iterator[Path]:
    """
    Materialize an AtomArray as a temporary file for file-based applications.

    Parameters
    ----------
    temp_file_format:
        Preferred alias for selecting temporary file format.
        If provided, it overrides file_format.

    The file is deleted after context exit by default.
    """
    resolved_format = _normalize_file_format(temp_file_format or file_format)
    resolved_suffix = suffix or f".{resolved_format}"
    if not resolved_suffix.startswith("."):
        resolved_suffix = f".{resolved_suffix}"

    fd, temp_name = tempfile.mkstemp(
        prefix=prefix,
        suffix=resolved_suffix,
        dir=Path(dir) if dir is not None else None,
    )
    os.close(fd)
    temp_path = Path(temp_name)

    try:
        write_structure_file(
            atom_array=atom_array,
            output_file=temp_path,
            file_format=resolved_format,
            **write_kwargs,
        )
        yield temp_path
    finally:
        if delete_file:
            temp_path.unlink(missing_ok=True)


def call_with_structure_file(
    atom_array: AtomArray,
    app_callable: Callable[..., T],
    *app_args,
    file_format: str = "pdb",
    temp_file_format: str | None = None,
    **app_kwargs,
) -> T:
    """
    Write AtomArray to a temporary file and pass its path to an app callable.

    The callable receives the temporary file path as the first positional argument.
    """
    resolved_format = temp_file_format or file_format
    with atom_array_as_temp_file(
        atom_array=atom_array,
        file_format=resolved_format,
    ) as temp_file:
        return app_callable(str(temp_file), *app_args, **app_kwargs)
