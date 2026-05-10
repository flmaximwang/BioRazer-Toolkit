import logging, sys, re


class NoDecorationFormatter(logging.Formatter):
    def format(self, record):
        record.msg = re.sub(r"\x1b\[[0-9;]*m", "", str(record.msg))
        if "\n" in record.msg:
            record.msg = "\n" + record.msg
        return super().format(record)


def get_default_console_logger(logger_name: str, level=logging.INFO) -> logging.Logger:
    logger = logging.getLogger(logger_name)
    # Check if the logger has handlers already to avoid duplicate logs
    formatter = NoDecorationFormatter(
        "[%(asctime)s @ %(name)s] %(levelname)s: %(message)s"
    )
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    console_handler.setLevel(level)
    logger.addHandler(console_handler)
    logger.setLevel(level)

    return logger


def get_default_file_logger(
    logger_name: str, log_file: str, level=logging.INFO
) -> logging.Logger:
    logger = logging.getLogger(logger_name)
    # Check if the logger has handlers already to avoid duplicate logs
    formatter = NoDecorationFormatter(
        "[%(asctime)s @ %(name)s] %(levelname)s: %(message)s"
    )
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    file_handler.setLevel(level)
    logger.addHandler(file_handler)
    logger.setLevel(level)

    return logger
