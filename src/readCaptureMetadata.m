function meta = readCaptureMetadata(metadataPath)
%readCaptureMetadata Read capture metadata from a JSON file.
%
%   meta = readCaptureMetadata(metadataPath) reads a JSON sidecar file and
%   returns a canonical metadata structure using the project field names.

arguments
    metadataPath (1,1) string
end

if ~isfile(metadataPath)
    error("readCaptureMetadata:FileNotFound", ...
        "Metadata file not found: %s", metadataPath);
end

rawText = fileread(metadataPath);
rawMeta = jsondecode(rawText);
[meta, ~] = validateMetadata(rawMeta);
meta.MetadataPath = char(metadataPath);

end
