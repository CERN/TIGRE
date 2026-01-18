"""This module is a trimmed down version of the image.py module from Pylinac for reading XIM images.
(source: https://github.com/jrkerns/pylinac/tree/master)."""

from __future__ import annotations

from pathlib import Path
from typing import BinaryIO
import numpy as np
import struct

XIM_PROP_INT = 0
XIM_PROP_DOUBLE = 1
XIM_PROP_STRING = 2
XIM_PROP_DOUBLE_ARRAY = 4
XIM_PROP_INT_ARRAY = 5


def decode_binary(
    file: BinaryIO,
    dtype: type[int] | type[float] | type[str] | str | np.dtype,
    num_values: int = 1,
    cursor_shift: int = 0,
    strip_empty: bool = True,
) -> int | float | str | np.ndarray | list:
    """Read in a raw binary file and convert it to given data types.

    Parameters
    ----------
    file
        The open file object.
    dtype
        The expected data type to return. If int or float and num_values > 1, will return numpy array.
    num_values
        The expected number of dtype to return

        .. note:: This is not the same as the number of bytes.

    cursor_shift : int
        The number of bytes to move the cursor forward after decoding. This is used if there is a
        reserved section after the read-in segment.
    strip_empty : bool
        Whether to strip trailing empty/null values for strings.
    """
    f = file

    if isinstance(dtype, str):
        s = struct.calcsize(dtype) * num_values
        output = struct.unpack(dtype * num_values, f.read(s))
        if len(output) == 1:
            output = output[0]
    elif dtype is str:
        ssize = struct.calcsize("c") * num_values
        output = struct.unpack("c" * num_values, f.read(ssize))
        if strip_empty:
            output = "".join(o.decode() for o in output if o != b"\x00")
        else:
            output = "".join(o.decode() for o in output)
    elif dtype is int:
        ssize = struct.calcsize("i") * num_values
        output = np.asarray(struct.unpack("i" * num_values, f.read(ssize)))
        if len(output) == 1:
            output = int(np.squeeze(output))
    elif dtype is float:
        ssize = struct.calcsize("f") * num_values
        output = np.asarray(struct.unpack("f" * num_values, f.read(ssize)))
        if len(output) == 1:
            output = float(np.squeeze(output))
    else:
        raise TypeError(f"datatype '{dtype}' was not valid")

    # shift cursor if need be (e.g. if a reserved section follows)
    if cursor_shift:
        f.seek(cursor_shift, 1)
    return output


class XIM:
    """A class to open, read, and/or export an .xim image, Varian's custom image format which is 99.999% PNG

    This had inspiration from a number of places:
    - https://gist.github.com/1328/7da697c71f9c4ef12e1e
    - https://medium.com/@duhroach/how-png-works-f1174e3cc7b7
    - https://www.mathworks.com/matlabcentral/answers/419228-how-to-write-for-loop-and-execute-data
    - https://www.w3.org/TR/PNG-Filters.html
    - https://bitbucket.org/dmoderesearchtools/ximreader/src/master/
    """

    array: np.ndarray
    properties: dict

    def __init__(self, file_path: str | Path, read_pixels: bool = True):
        """
        Parameters
        ----------
        file_path
            The path to the file of interest.
        read_pixels
            Whether to read and parse the pixel information. Doing so is quite slow.
            Set this to false if, e.g., you are searching for images only via tags or doing
            a pre-filtering of image selection.
        """

        with open(file_path, "rb") as xim:
            self.format_id = decode_binary(xim, str, 8)
            self.format_version = decode_binary(xim, int)
            self.img_width_px = decode_binary(xim, int)
            self.img_height_px = decode_binary(xim, int)
            self.bits_per_pixel = decode_binary(xim, int)
            self.bytes_per_pixel = decode_binary(xim, int)
            self.compression = decode_binary(xim, int)
            if not self.compression:
                pixel_buffer_size = decode_binary(xim, int)
                self.pixel_buffer = decode_binary(xim, str, num_values=pixel_buffer_size)
            else:
                lookup_table_size = decode_binary(xim, int)
                self.lookup_table = np.fromfile(xim, count=lookup_table_size, dtype=np.uint8)
                if read_pixels:
                    lookup_keys = self._parse_lookup_table(self.lookup_table)
                    self.array = self._parse_compressed_bytes(xim, lookup_table=lookup_keys)
                else:
                    comp_pixel_buffer_size = decode_binary(xim, int)
                    _ = decode_binary(xim, "c", num_values=comp_pixel_buffer_size)
                decode_binary(xim, int)
            self.num_hist_bins = decode_binary(xim, int)
            self.histogram = decode_binary(xim, int, num_values=self.num_hist_bins)
            self.num_properties = decode_binary(xim, int)
            self.properties = {}
            for prop in range(self.num_properties):
                name_length = decode_binary(xim, int)
                name = decode_binary(xim, str, num_values=name_length)
                tipe = decode_binary(xim, int)
                if tipe == XIM_PROP_INT:
                    value = decode_binary(xim, int)
                elif tipe == XIM_PROP_DOUBLE:
                    value = decode_binary(xim, "d")
                elif tipe == XIM_PROP_STRING:
                    num_bytes = decode_binary(xim, int)
                    value = decode_binary(xim, str, num_values=num_bytes)
                elif tipe == XIM_PROP_DOUBLE_ARRAY:
                    num_bytes = decode_binary(xim, int)
                    value = decode_binary(
                        xim, "d", num_values=int(num_bytes // 8)
                    )  # doubles are 8 bytes
                elif tipe == XIM_PROP_INT_ARRAY:
                    num_bytes = decode_binary(xim, int)
                    value = decode_binary(
                        xim, int, num_values=int(num_bytes // 4)
                    )  # ints are 4 bytes
                self.properties[name] = value

    @staticmethod
    def _parse_lookup_table(lookup_table_bytes: np.ndarray) -> np.ndarray:
        """The lookup table doesn't follow normal structure conventions like 1, 2, or 4 byte values. They
        got smart and said each value is 2 bits. Yes, bits. This means each byte is actually 4 values.
        Python only reads things as granular as bytes. To get around this the general logic is:

        1) interpret the data as integers at the single byte level
        2) convert those integers back into bit representation; e.g. 115 => 01110011. Note the representation must contain the full byte. I.e. 3 => 11 does not work.
        3) split the binary representation into the 2-bit representations; generates 4x the number of elements. 01110011 => (01, 11, 00, 11)
        4) Convert the 2-bit representation back into integers (01, 11, 00, 11) => (1, 3, 0, 3)

        .. note::

            This is ripe for optimization, but brevity and clarity won out. Options include bit-shifting (fastest)
            and numpy.packbits/unpackbits.
        """
        bit_shift = np.array([0, 2, 4, 6])
        lookup_table = (lookup_table_bytes[:, np.newaxis] >> bit_shift[np.newaxis, :]) & 0b00000011
        return lookup_table.flatten()

    def _parse_compressed_bytes(self, xim: BinaryIO, lookup_table: np.ndarray) -> np.ndarray:
        """Parse the compressed pixels. We have to do this pixel-by-pixel because each
        pixel can have a different number of bytes representing it

        Per the readme:

        1) The first row is uncompressed
        2) The first element of the second row is uncompressed
        3) all other elements are represented by 1, 2, or 4 bytes of data (the annoying part)
        4) The byte size of the element is given in the lookup table

        So, we have to read in 1, 2, or 4 bytes and convert to an integer depending on
        the lookup table, which tells us how many bytes to read in

        .. note::

            Optimization can help here. A few ideas:

            - reading in groups of data of the same byte size. I already tried this, and I think it will work, but I couldn't get it going.
            - reading in rows of data where no byte change occurred in that row. Similar to above.
            - Using joblib or a processpool
        """
        img_height = self.img_height_px
        img_width = self.img_width_px
        if self.bytes_per_pixel == 1:
            dtype = np.int8
        elif self.bytes_per_pixel == 2:
            dtype = np.int16
        elif self.bytes_per_pixel == 4:
            dtype = np.int32
        elif self.bytes_per_pixel == 8:
            dtype = np.int64
        else:
            raise ValueError(
                "The XIM image has an unsupported bytes per pixel value. Raise a ticket on the pylinac Github with this file."
            )
        # first row and 1st element, 2nd row is uncompressed
        # this SHOULD work by reading the # of bytes specified in the header but AFAICT this is just a standard int (4 bytes)
        compressed_array = self._get_diffs(lookup_table, xim, dtype, img_width, img_height)

        compressed_array = compressed_array.reshape((img_height, img_width))
        compressed_array_windows = np.lib.stride_tricks.sliding_window_view(
            compressed_array, (2, img_width), writeable=True
        )

        next_row_correction = 0
        for window in compressed_array_windows[:, 0]:
            window[1, 0] = np.add(window[1, 0], next_row_correction)
            window[1, 1:] += window[0, 1:] - window[0, :-1]
            window[1] = np.cumsum(window[1])

            next_row_correction = np.subtract(np.add(window[1, -1], window[1, 0]), window[0, -1])

        return compressed_array

    @staticmethod
    def _get_diffs(
        lookup_table: np.ndarray,
        xim: BinaryIO,
        dtype: np.dtype,
        img_width: int,
        img_height: int,
    ):
        """Read in all the pixel value 'diffs'. These can be 1, 2, or 4 bytes in size,
        so instead of just reading N pixels of M bytes which would be SOOOO easy, we have to read dynamically

        We optimize here by reading bytes in clumps, which is way faster than reading one at a time.
        Knowing that most values are single bytes with an occasional 2-byte element
        we read chunks that all look like (n 1-bytes and 1 2-byte)
        """
        comp_pixel_buffer_size = decode_binary(xim, int)
        file_array = np.fromfile(xim, dtype=np.uint8, count=comp_pixel_buffer_size)

        compressed_array = np.zeros((img_height * img_width), dtype=dtype)
        compressed_array[: img_width + 1] = file_array[: (img_width + 1) * 4].view(np.int32)
        file_array = file_array[(img_width + 1) * 4 :]

        change_indices = np.where(np.diff(lookup_table) != 0)[0] + 1
        lengths = np.diff(np.concatenate(([0], change_indices, [len(lookup_table)])))
        values = lookup_table[np.concatenate(([0], change_indices))]

        len_diffs = img_width * img_height - img_width - 1
        LOOKUP_CONVERSION = {0: "<i1", 1: "<i2", 2: "<i4"}
        start = 0
        for value, length in zip(values, lengths):
            read_dtype = LOOKUP_CONVERSION[value]
            length = min(length, len_diffs - start)
            stop = start + length
            bytes_len = length * (1 << value)
            compressed_array[img_width + 1 + start : img_width + 1 + stop] = file_array[
                :bytes_len
            ].view(read_dtype)
            file_array = file_array[bytes_len:]
            start = stop
        return compressed_array
