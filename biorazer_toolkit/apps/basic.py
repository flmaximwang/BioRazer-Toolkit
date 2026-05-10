from abc import abstractmethod
from pathlib import Path
from dataclasses import dataclass
import codecs, errno, logging, subprocess, selectors, pty, os, sys, re


@dataclass
class App:
    """
    Base class for all applications. Provides common methods for running subprocesses and logging.

        Attributes
        ----------
        dir : Path
            The directory where the application is installed.
        bin : Path
            The path to the application's executable. Rosetta has multiple executables, so this can be set to any of them to use different tools.
        logger : logging.Logger
            Logger for the application.
            - If not provided, a default logger will be created.
            - If a string is provided, it will be used as the logger name to create a logger.
            - You can call set_default_logger_style() to set the default logging style, level and handlers.
    """

    app_dir: str | Path = None
    app_bin: str | Path = None
    logger: logging.Logger = None

    def __post_init__(self):
        if isinstance(self.logger, str):
            self.logger = logging.getLogger(self.logger)
        if self.logger is None:
            self.logger = logging.getLogger(__name__)
            self.set_default_logger_style()
        self.app_dir = Path(self.app_dir) if self.app_dir else None
        self.app_bin = Path(self.app_bin) if self.app_bin else None
        self.logger.debug(
            f"Post-init App with dir: {self.app_dir}, bin: {self.app_bin}, logger: {self.logger}"
        )

    def set_default_handler(
        self,
        handler_types=("StreamHandler",),
        file_path=None,
        stream=None,
        file_mode="a",
        file_encoding="utf-8",
    ):
        """
        Parameters
        ----------
        - handler_types: str or list of str or logging.Handler subclass or list of logging.Handler subclasses
            Types of handlers to set as default. Supported string values are "StreamHandler" and "FileHandler".
        - stream: file-like object
            Stream to use for StreamHandler. Defaults to sys.stdout if not provided.
        """
        # Remove existing handlers before applying new defaults.
        for handler in list(self.logger.handlers):
            self.logger.removeHandler(handler)
            handler.close()

        if isinstance(handler_types, (str, type)):
            handler_types = [handler_types]

        for handler_type in handler_types:
            if isinstance(handler_type, str):
                handler_key = handler_type.lower()
                if handler_key in {"streamhandler", "stream"}:
                    handler = logging.StreamHandler(
                        stream if stream is not None else sys.stdout
                    )
                elif handler_key in {"filehandler", "file"}:
                    if file_path is None:
                        raise ValueError("file_path is required for FileHandler.")
                    handler = logging.FileHandler(
                        filename=str(file_path),
                        mode=file_mode,
                        encoding=file_encoding,
                    )
                else:
                    raise ValueError(f"Unsupported handler type: {handler_type}")
            elif isinstance(handler_type, type) and issubclass(
                handler_type, logging.FileHandler
            ):
                if file_path is None:
                    raise ValueError("file_path is required for FileHandler.")
                handler = handler_type(
                    filename=str(file_path),
                    mode=file_mode,
                    encoding=file_encoding,
                )
            elif isinstance(handler_type, type) and issubclass(
                handler_type, logging.StreamHandler
            ):
                handler = handler_type(stream if stream is not None else sys.stdout)
            else:
                raise ValueError(f"Unsupported handler type: {handler_type}")

            self.logger.addHandler(handler)

    def set_default_logger_style(
        self,
        fmt="[%(asctime)s @ %(name)s] %(levelname)s: %(message)s",
        level=logging.INFO,
        handler_types=("StreamHandler",),
        file_path=None,
    ):
        formatter = logging.Formatter(fmt)

        self.set_default_handler(handler_types=handler_types, file_path=file_path)
        for handler in self.logger.handlers:
            handler.setFormatter(formatter)
            handler.setLevel(level)
        self.logger.setLevel(level)

    def _normalize_log_message(self, message):
        if message is None:
            return ""
        if not isinstance(message, str):
            message = str(message)

        # Strip ANSI escape sequences used for terminal styling.
        return re.sub(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])", "", message)

    @abstractmethod
    def from_default_bin(cls, bin_name):
        app = cls(app_dir=None, app_bin=bin_name)
        return app

    @abstractmethod
    def set_dir(self, app_dir: str | Path):
        self.app_dir = Path(app_dir)

    @abstractmethod
    def set_bin(self, app_bin: str | Path):
        self.app_bin = Path(app_bin)

    @abstractmethod
    def check(self):
        """
        Implement in subclass to check if the app is installed correctly.
        """

    @abstractmethod
    def run_with_structure(
        self,
        atom_array,
        *args,
        input_file_format="pdb",
        **kwargs,
    ):
        """
        Run app from an in-memory structure by materializing a temporary file.

        Subclasses should implement app-specific argument wiring.
        """
        raise NotImplementedError("Subclasses must implement run_with_structure().")

    def run(self, *args, cwd=".", get_output=True, verbose=True, mode="pty"):
        """
        Run self.bin with given args and kwargs.


        Parameters
        ----------
        args : list
            List of command line arguments to pass to the application. "-s" and the input file path should not be included in one "" argument, but should be passed as separate arguments.
        mode : str
            Mode to run the subprocess.
            - pty: Pseudo-terminal mode to ensure line-buffered output.
            - subprocess.run: Using subprocess.run to capture output.
            - subprocess.Popen: Using subprocess.Popen to capture output.
        """
        cmd_args = [f"{self.app_bin}"]
        for arg in args:
            cmd_args.append(str(arg))
        self.logger.info(f"Running command: {' '.join(cmd_args)}")

        if mode == "pty":
            # 通过 PTY 确保 C++ 的输出采用行缓冲, 而不是全缓冲
            master_stdout, slave_stdout = pty.openpty()
            master_stderr, slave_stderr = pty.openpty()

            p = subprocess.Popen(
                cmd_args,
                cwd=cwd,
                stdout=slave_stdout,
                stderr=slave_stderr,
                close_fds=True,
            )
            os.close(slave_stdout)  # Prevent child from hanging on to the fd
            os.close(slave_stderr)

            output = ""
            error = ""
            finish = False
            sel = selectors.DefaultSelector()
            for master_fd in (master_stdout, master_stderr):
                sel.register(
                    master_fd, selectors.EVENT_READ
                )  # Register for read events

            line_left = {
                master_stdout: "",
                master_stderr: "",
            }
            decoders = {
                master_stdout: codecs.getincrementaldecoder("utf-8")(),
                master_stderr: codecs.getincrementaldecoder("utf-8")(),
            }
            while True:
                self.logger.debug("Process is still running...")

                self.logger.debug(f"Checking if process has finished: {p.poll()}")
                if p.poll() is not None:
                    finish = True

                events = sel.select(timeout=1.0)  # Wait for events with a timeout
                self.logger.debug(f"Events received: {events}")

                if finish and len(events) == 0:
                    break  # Process finished and no more events

                for key, mask in events:
                    if (
                        mask & selectors.EVENT_READ
                    ):  # Bitwise AND to check whether mask contains EVENT_READ
                        try:
                            data = os.read(key.fileobj, 8192)
                        except OSError as exc:
                            # On Linux, reading from a PTY master after the slave
                            # has closed can raise EIO instead of returning b"".
                            if exc.errno == errno.EIO:
                                data = b""
                            else:
                                raise
                        if not data:
                            remaining = decoders[key.fileobj].decode(b"", final=True)
                            remaining = line_left[key.fileobj] + remaining
                            line_left[key.fileobj] = ""
                            if remaining:
                                normalized_line = self._normalize_log_message(
                                    remaining
                                ).strip()
                                if key.fileobj == master_stderr:
                                    if verbose:
                                        self.logger.error(normalized_line)
                                    if get_output:
                                        output += remaining
                                    error += remaining
                                elif key.fileobj == master_stdout:
                                    if verbose:
                                        self.logger.info(normalized_line)
                                    if get_output:
                                        output += remaining
                            sel.unregister(key.fileobj)  # No more data to read
                            continue

                        decoded = decoders[key.fileobj].decode(data, final=False)
                        if not decoded:
                            continue
                        lines_tmp = decoded.splitlines(keepends=True)
                        self.logger.debug(f"Line left: {line_left[key.fileobj]}")
                        self.logger.debug(f"Read data: {lines_tmp}")
                        if line_left[key.fileobj]:
                            lines_tmp[0] = line_left[key.fileobj] + lines_tmp[0]
                            line_left[key.fileobj] = ""
                        if lines_tmp[-1].endswith("\n"):
                            line_left[key.fileobj] = ""
                        else:
                            line_left[key.fileobj] = lines_tmp.pop()
                        self.logger.debug(f"Line left: {line_left[key.fileobj]}")
                        for line in lines_tmp:
                            normalized_line = self._normalize_log_message(line).strip()
                            if key.fileobj == master_stderr:
                                if verbose:
                                    self.logger.error(normalized_line)
                                if get_output:
                                    output += line
                                error += line
                            elif key.fileobj == master_stdout:  # master_stdout
                                if verbose:
                                    self.logger.info(normalized_line)
                                if get_output:
                                    output += line
                            else:
                                raise RuntimeError(f"Unknown file object {key.fileobj}")
                    else:
                        raise RuntimeError(f"Unknown event mask {mask}")

            if p.returncode != 0:
                if error:
                    raise RuntimeError(f"{self.app_bin} failed with error:\n{error}")
                else:
                    raise RuntimeError(f"{self.app_bin} failed with output\n{output}")

            os.close(master_stdout)
            os.close(master_stderr)
            if get_output:
                return output
        elif mode == "subprocess.run":
            result = subprocess.run(
                cmd_args,
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode != 0:
                if result.stderr:
                    raise RuntimeError(
                        f"{self.app_bin} failed with error:\n{result.stderr}"
                    )
                else:
                    raise RuntimeError(
                        f"{self.app_bin} failed with output\n{result.stdout}"
                    )

            if result.stdout:
                output = result.stdout
            else:
                output = result.stderr
            if verbose:
                self.logger.info("\n" + self._normalize_log_message(output))
            if get_output:
                return output
        elif mode == "subprocess.Popen":
            p = subprocess.Popen(
                cmd_args,
                bufsize=1,
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            stdout = ""
            stderr = ""
            while p.poll() is None:
                self.logger.debug("Process is still running...")
                for line in p.stdout:
                    if verbose:
                        self.logger.info(self._normalize_log_message(line).strip())
                    stdout += line
                for line in p.stderr:
                    if verbose:
                        self.logger.error(self._normalize_log_message(line).strip())
                    stderr += line

            stderr += p.stderr.read()
            stdout += p.stdout.read()

            if p.returncode != 0:
                if stderr:
                    raise RuntimeError(f"{self.app_bin} failed with error:\n{stderr}")
                else:
                    raise RuntimeError(f"{self.app_bin} failed with output\n{stdout}")

            if stdout:
                output = stdout
            else:
                output = stderr
            if verbose:
                self.logger.info("\n" + self._normalize_log_message(output))
            if get_output:
                return output

        else:
            raise ValueError(f"Unknown mode: {mode}")
