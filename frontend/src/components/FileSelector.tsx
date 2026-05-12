import { open } from "@tauri-apps/plugin-dialog";
import { useState } from "react";

/** File-type icons — colors based on biological meaning */
function FileIcon({ type }: { type: string }) {
  const colors: Record<string, string> = {
    VCF: "#E52421",        /* Red — variants / mutations */
    Variants: "#E52421",   /* Red — variants / mutations */
    FASTA: "#3460AA",      /* Blue — reference sequence */
    Annotation: "#66B32E", /* Green — genes */
    BAM: "#EF781A",        /* Orange — alignments */
  };
  const color = colors[type] ?? "#888";
  return (
    <svg width="38" height="38" viewBox="0 0 32 32" fill="none">
      <path
        d="M8 4h10l8 8v16a2 2 0 01-2 2H8a2 2 0 01-2-2V6a2 2 0 012-2z"
        fill={color}
        opacity="0.12"
        stroke={color}
        strokeWidth="1.5"
      />
      <path d="M18 4v8h8" stroke={color} strokeWidth="1.5" />
      <text
        x="14"
        y="24"
        textAnchor="middle"
        fill={color}
        fontSize="6"
        fontWeight="700"
        fontFamily="system-ui"
      >
        {type.slice(0, 4)}
      </text>
    </svg>
  );
}

function CheckIcon() {
  return (
    <svg width="18" height="18" viewBox="0 0 18 18" fill="none">
      <circle cx="9" cy="9" r="9" fill="#149389" opacity="0.15" />
      <path
        d="M5.5 9.5l2 2 5-5"
        stroke="#149389"
        strokeWidth="1.8"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
    </svg>
  );
}

function ClearIcon() {
  return (
    <svg width="16" height="16" viewBox="0 0 16 16" fill="none">
      <path
        d="M4 4l8 8M12 4l-8 8"
        stroke="currentColor"
        strokeWidth="1.5"
        strokeLinecap="round"
      />
    </svg>
  );
}

function UploadIcon() {
  return (
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none">
      <path
        d="M12 16V4m0 0l-4 4m4-4l4 4"
        stroke="currentColor"
        strokeWidth="1.8"
        strokeLinecap="round"
        strokeLinejoin="round"
      />
      <path
        d="M4 17v2a2 2 0 002 2h12a2 2 0 002-2v-2"
        stroke="currentColor"
        strokeWidth="1.8"
        strokeLinecap="round"
      />
    </svg>
  );
}

interface FileSelectorProps {
  label: string;
  value: string;
  onChange: (path: string) => void;
  filters?: { name: string; extensions: string[] }[];
  required?: boolean;
  isDragActive?: boolean;
  justAssigned?: boolean;
  fileType?: string;
  validationError?: boolean;
}

/** Extract filename from a full path */
function basename(path: string): string {
  return path.split(/[\\/]/).pop() ?? path;
}

export default function FileSelector({
  label,
  value,
  onChange,
  filters,
  required,
  isDragActive,
  justAssigned,
  fileType,
  validationError,
}: FileSelectorProps) {
  const [hovering, setHovering] = useState(false);
  const hasFile = value.length > 0;
  const filterName = filters?.[0]?.name ?? "File";

  const handleBrowse = async () => {
    const selected = await open({
      multiple: false,
      filters: filters,
    });
    if (selected) {
      onChange(selected as string);
    }
  };

  const handleClear = (e: React.MouseEvent) => {
    e.stopPropagation();
    onChange("");
  };

  const zoneClass = [
    "drop-zone",
    hasFile ? "drop-zone--filled" : "",
    isDragActive ? "drop-zone--drag-active" : "",
    hovering ? "drop-zone--hover" : "",
    justAssigned ? "drop-zone--just-assigned" : "",
    validationError && !hasFile ? "drop-zone--error" : "",
  ]
    .filter(Boolean)
    .join(" ");

  return (
    <div className="file-selector">
      <div className="file-selector-label">
        <span className="file-selector-label-text">
          {label}
          {required && <span className="required">*</span>}
        </span>
        <span className="file-selector-extensions">
          {filters?.[0]?.extensions.map((e) => `.${e}`).join(", ")}
        </span>
      </div>
      <div
        className={zoneClass}
        data-type={fileType ?? filterName}
        onClick={handleBrowse}
        onMouseEnter={() => setHovering(true)}
        onMouseLeave={() => setHovering(false)}
        role="button"
        tabIndex={0}
        onKeyDown={(e) => e.key === "Enter" && handleBrowse()}
      >
        {hasFile ? (
          <div className="drop-zone-filled">
            <FileIcon type={filterName} />
            <div className="drop-zone-file-info">
              <div className="drop-zone-filename">
                <CheckIcon />
                <span>{basename(value)}</span>
              </div>
              <div className="drop-zone-path" title={value}>
                {value}
              </div>
            </div>
            <button
              className="drop-zone-clear"
              onClick={handleClear}
              title="Remove file"
              type="button"
            >
              <ClearIcon />
            </button>
          </div>
        ) : (
          <div className="drop-zone-empty">
            <UploadIcon />
            <div className="drop-zone-empty-text">
              <span className="drop-zone-empty-main">
                Drop <strong>{filterName}</strong> file here
              </span>
              <span className="drop-zone-empty-sub">
                or click to browse
              </span>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
