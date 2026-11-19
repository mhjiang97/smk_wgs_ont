#!/usr/bin/env python3
"""Run SVision 1.3.8 without forking a TensorFlow/OpenCV prediction worker.

Collection retains SVision's worker pool. Prediction runs in the main process,
with TensorFlow imported after collection. The model and calling parameters
are unchanged. --resume-segments is for explicitly verified, complete segment
files from an interrupted run of the same sample and parameters.
"""
import ast
import os
from pathlib import Path


def assignment(node, name):
    return (isinstance(node, ast.Assign) and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name) and node.targets[0].id == name)


def find_assignment(body, name):
    matches = [i for i, node in enumerate(body) if assignment(node, name)]
    if len(matches) != 1:
        raise RuntimeError('Unsupported SVision source: expected one assignment to ' + name)
    return matches[0]


def adapted_tree(source):
    tree = ast.parse(source)
    main_blocks = [node for node in tree.body if isinstance(node, ast.If)
                   and isinstance(node.test, ast.Compare)
                   and isinstance(node.test.left, ast.Name)
                   and node.test.left.id == '__name__']
    if len(main_blocks) != 1:
        raise RuntimeError('Unsupported SVision main block')
    main = main_blocks[0]
    imports = [node for node in tree.body if isinstance(node, ast.ImportFrom)
               and node.module == 'src.network.predict']
    predictors = [node for node in main.body if isinstance(node, ast.FunctionDef)
                  and node.name == 'predict_one_chrom']
    if len(imports) != 1 or len(predictors) != 1:
        raise RuntimeError('Unsupported SVision prediction implementation')
    # Import TensorFlow only after the collection workers have exited.
    tree.body.remove(imports[0])
    predictors[0].body.insert(0, imports[0])

    pools = [i for i, node in enumerate(main.body)
             if assignment(node, 'process_pool') and isinstance(node.value, ast.Call)
             and isinstance(node.value.func, ast.Attribute)
             and isinstance(node.value.func.value, ast.Name)
             and node.value.func.value.id == 'multiprocessing'
             and node.value.func.attr == 'Pool']
    if len(pools) != 2:
        raise RuntimeError('Unsupported SVision worker pools')
    prediction_end = find_assignment(main.body, 'end_time_predict')
    if not pools[1] < prediction_end:
        raise RuntimeError('Unsupported SVision prediction block')
    main.body[pools[1]:prediction_end] = ast.parse('''
for chrom in task_list_bychrom.keys():
    logging.info('[Prediction process] Running in the main process without fork')
    prediction_error = predict_one_chrom(chrom, predict_results_dir, options)
    if prediction_error is not None:
        raise RuntimeError(prediction_error)
''').body

    parsers = [node for node in tree.body if isinstance(node, ast.FunctionDef)
               and node.name == 'parse_arguments']
    if len(parsers) != 1:
        raise RuntimeError('Unsupported SVision argument parser')
    parser = parsers[0]
    options_at = find_assignment(parser.body, 'options')
    parser.body[options_at:options_at] = ast.parse('''
parser.add_argument('--resume-segments', action='store_true',
    help='Reuse verified complete segments for a single chromosome')
''').body

    collection_start = find_assignment(main.body, 'chrom_split_files')
    collection_end = find_assignment(main.body, 'end_time_collect')
    collection = main.body[collection_start:collection_end]
    resume = ast.parse('''
if options.resume_segments:
    if len(task_list_bychrom) != 1 or options.chrom not in task_list_bychrom:
        raise RuntimeError('--resume-segments requires a single complete chromosome')
    for chrom in task_list_bychrom:
        segment_file = os.path.join(segments_out_path, chrom + '.segments.all.bed')
        if not os.path.isfile(segment_file):
            raise RuntimeError('Missing completed segment file: ' + segment_file)
        logging.info('[Resume] Using verified completed segments: ' + segment_file)
else:
    pass
''').body[0]
    resume.orelse = collection
    main.body[collection_start:collection_end] = [resume]
    return ast.fix_missing_locations(tree)


def main():
    source = Path('/usr/local/bin/SVision')
    tree = adapted_tree(source.read_text())
    os.environ['OPENCV_OPENCL_RUNTIME'] = 'disabled'
    import cv2
    cv2.ocl.setUseOpenCL(False)
    cv2.setNumThreads(1)
    namespace = {'__name__': '__main__', '__file__': str(source)}
    exec(compile(tree, str(source), 'exec'), namespace)


if __name__ == '__main__':
    main()
