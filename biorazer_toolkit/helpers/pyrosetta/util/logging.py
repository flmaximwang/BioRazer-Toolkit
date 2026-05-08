import pyrosetta
import logging, re


class PyRosettaFormatter(logging.Formatter):

    def format(self, record: logging.LogRecord) -> str:
        record.msg = re.sub(r"\x1b\[[0-9;]*m", "", record.msg)
        if "\n" in record.msg:
            record.msg = "\n" + record.msg
        else:
            record.msg = record.msg + " ;"
        return super().format(record)


def setup_logging(log_file: str, level=logging.INFO, mode="a") -> None:
    """
    Run this command before pyrosetta.init() to set up logging to a file.
    This will log all PyRosetta output to the specified log_file.
    """

    logger = logging.getLogger("rosetta")
    logger.setLevel(level)

    # Create file handler which logs even debug messages
    fh = logging.FileHandler(log_file, mode=mode)
    fh.setLevel(level)

    # Create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)

    # Create formatter and add it to the handlers
    formatter = PyRosettaFormatter(
        "[ %(asctime)s @ %(name)s ] %(levelname)s: %(message)s"
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # Add the handlers to the logger
    if logger.hasHandlers():
        logger.handlers.clear()
    logger.addHandler(fh)
    logger.addHandler(ch)

    pyrosetta.logging_support.set_logging_sink()
