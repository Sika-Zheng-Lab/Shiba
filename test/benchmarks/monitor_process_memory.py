"""Linux process-tree memory sampler used in the 2026-10-15 validation.

Optional dependency: psutil. Outputs pipeline.log, memory.tsv and
memory_summary.json; use a new directory for each measurement.
Phase detection is specific to Shiba beta_binomial_glm.R command lines.
"""
import argparse, json, os, subprocess, time
from pathlib import Path
import psutil
p = argparse.ArgumentParser()
p.add_argument('--directory', required=True)
p.add_argument('command', nargs=argparse.REMAINDER)
a = p.parse_args()
directory = Path(a.directory)
directory.mkdir(parents=True, exist_ok=True)
if not a.command:
    p.error('a command is required after --')
command = a.command[1:] if a.command[0] == '--' else a.command
if not command:
    p.error('a command is required after --')
start = time.time()
peaks = {}
fit_seen = False
errors = 0
samples = 0
max_n = 0
with (directory / 'pipeline.log').open('w') as logfile, (directory / 'memory.tsv').open('w') as memory:
    proc = subprocess.Popen(command, stdout=logfile, stderr=subprocess.STDOUT)
    memory.write('seconds\tphase\tprocesses\tpss_bytes\tuss_bytes\trss_bytes\tpython_pss_bytes\tr_pss_bytes\n')
    while proc.poll() is None:
        try:
            children = [psutil.Process(proc.pid)] + psutil.Process(proc.pid).children(recursive=True)
        except psutil.NoSuchProcess:
            break
        total = dict(pss=0, uss=0, rss=0, python_pss=0, r_pss=0)
        r_stage = None
        n = 0
        for child in children:
            try:
                info = child.memory_full_info()
                cmd = child.cmdline()
                n += 1
                for key in ('pss', 'uss', 'rss'):
                    total[key] += getattr(info, key)
                is_r = any(('beta_binomial_glm.R' in arg for arg in cmd))
                total['r_pss' if is_r else 'python_pss'] += info.pss
                if is_r:
                    r_stage = 'fitting' if 'fit' in cmd else 'design'
            except psutil.NoSuchProcess:
                pass
            except psutil.AccessDenied:
                errors += 1
        if r_stage == 'fitting':
            fit_seen = True
        phase = r_stage or ('post_fit' if fit_seen else 'input_and_preparation')
        for label in ('overall', phase):
            peak = peaks.setdefault(label, {k: 0 for k in total})
            for key in total:
                peak[key] = max(peak[key], total[key])
        memory.write(f'{time.time() - start:.3f}\t{phase}\t{n}\t' + '\t'.join((str(total[k]) for k in total)) + '\n')
        memory.flush()
        samples += 1
        max_n = max(max_n, n)
        time.sleep(1)
    code = proc.wait()
report = {'command': command, 'exit_code': code, 'elapsed_seconds': time.time() - start, 'peak_bytes': peaks, 'samples': samples, 'memory_access_errors': errors, 'maximum_processes': max_n, 'cpu_affinity': sorted(os.sched_getaffinity(0)), 'slurm_job_id': os.environ.get('SLURM_JOB_ID'), 'host': os.uname().nodename, 'sampling': 'One second sleep plus process-tree memory scan; sum PSS avoids shared-page double counting.'}
(directory / 'memory_summary.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report), flush=True)
raise SystemExit(code)
