#!/bin/bash
# 批量下载ALIEN存储中的AnalysisResults.root文件（边检测边下载版）

# ===================== HY编号列表 =====================
HY_NUMBERS=(
5109428 5109427 5109426 5109425
5109424 5109423 5109422 5109421
5109420 5109419 5109418 5109417
5109416 5109415 5109414 5109413
5109412 5109411 5109410 5109409
5109408 5109407 5109406 5109405
5109404 5109403 5109402 5109401
5109400 5109399 5109398 5109397
5109396 5109395 5109394 5109393
5109392 5109391 5109390 5109389
5109388 5109387 5109386 5109385
5109384 5109383 5109382 5109381
5109380 5109379 5109378 5109377
5109376 5109375 5109374 5109373
5109372 5109371 5109370 5109369
5109368 5109367 5109366 5109365
5109364 5109363 5109362 5109361
5109360 5109359 5109358 5109357
5109356 5109355 5109354 5109353
5109350 5109349 5109348 5109347
5109346 5109345 5109344 5109343
5109342 5109341 5109340 5109339
5109338 5109337 5109336 5109335
5109334 5109333 5109332 5109331
5109326 5109325 5109324 5109323
5109322 5109321 5109320 5109319
5109318 5109317 5109316 5109315
5109314 5109313 5109312 5109311
5109310 5109309 5109308 5109306
5109305 5109304 5109303 5109302
5109301 5109300 5109299 5109298
5109297 5109296 5109295 5109294
5109293 5109292 5109291 5109290
5109289 5109288 5109287 5109286
5109285 5109284 5109283 5109282
5109281 5109280 5109279 5109278
5109277 5109276 5109275 5109274
5109273 5109272 5109271 5109270
5109269 5109268 5109267 5109266
5109265 5109264 5109263 5109262
5109261 5109260 5109259 5109258
5109257 5109256 5109255 5109254
5109253 5109252 5109251 5109250
5109249 5109248 5109247 5109246
5109245 5109244 5109243 5109242
5109241 5109240 5109239 5109238
5109237 5109236 5109235 5109234
5109233 5109232 5109231 5109230
5109229 5109228 5109227 5109226
5109225 5109224 5109223 5109222
5109221 5109220 5109219 5109218
5109217 5109216 5109215 5109214
5109213 5109212 5109211 5109210
5109209 5109208 5109207 5109206
5109205 5109204 5109203 5109202
5109201 5109200 5109199 5109198
5109197 5109196 5109195 5109194
5109193 5109192 5109191 5109190
5109189 5109188 5109187 5109186
5109185 5109184 5109183 5109182
5109181 5109180 5109179 5109178
5109177 5109176 5109175 5109174
5109173
)


# ===================== 核心配置 =====================
ALIEN_BASE_PATH="/alice/cern.ch/user/a/alihyperloop/jobs/0510"
LOCAL_SAVE_DIR="/home/suyoupeng/Downloads/lambdaspincorefile/Analysiscode/LHC23_pass4_thin_690954"
TIMEOUT=7200          # 单次下载超时（2小时）
SINGLE_FILE_RETRY=5   # 单个文件最多重试次数
THREADS=4             # 线程数（稳定优先）
LOG_FILE="${LOCAL_SAVE_DIR}/batch_download.log"
CHECK_INTERVAL=1      # 进度刷新间隔
DELAY_CHECK=10        # 下载后延迟判断时间

# ===================== 进度条配置 =====================
PROGRESS_BAR_LENGTH=50
GLOBAL_PROGRESS_LENGTH=60

# ===================== 全局统计变量 =====================
TOTAL_FILES=0          # 总文件数（动态累加）
DOWNLOADED_FILES=0
FAILED_FILES=0
SKIPPED_FOLDERS=0
RETRY_COUNT=0

# ===================== 检查alien环境 =====================
check_alien_env() {
    if ! command -v alien_ls &> /dev/null; then
        echo "❌ alien_ls 命令未找到，请执行: alienv enter AliPhysics/latest"
        exit 1
    fi
    if ! command -v alien_cp &> /dev/null; then
        echo "❌ alien_cp 命令未找到，请执行: alienv enter AliPhysics/latest"
        exit 1
    fi
    echo "✅ ALIEN环境检测通过"
}

# ===================== 检查ALIEN路径是否存在 =====================
check_alien_path_exists() {
    local alien_path=$1
    alien_ls "${alien_path}" > /dev/null 2>&1
    return $?
}

# ===================== 获取子文件夹列表 =====================
get_subfolders() {
    local alien_path="$1"
    
    alien_ls "${alien_path}" 2>/dev/null | while IFS= read -r line; do
        [[ -z "$line" ]] && continue
        
        # 检查是否以/结尾（表示文件夹）
        if [[ "$line" == */ ]]; then
            folder_name="${line%/}"
            if [[ "$folder_name" =~ ^[0-9]+$ ]]; then
                echo "$folder_name"
            fi
        fi
    done | sort
}

# ===================== 获取文件列表 =====================
get_files_in_folder() {
    local alien_path="$1"
    
    alien_ls "${alien_path}" 2>/dev/null | while IFS= read -r line; do
        [[ -z "$line" || "$line" == */ ]] && continue
        echo "$line"
    done
}

# ===================== 检查文件是否存在 =====================
check_alien_file_exists() {
    local alien_path=$1
    local dir_path=$(dirname "$alien_path")
    local file_name=$(basename "$alien_path")
    
    get_files_in_folder "$dir_path" | grep -qx "$file_name"
    return $?
}

# ===================== 进度条函数 =====================
show_progress() {
    local current=$1
    local total=$2
    local filename=$3

    local current_human=$(numfmt --to=iec --suffix=B --format="%.1f" $current 2>/dev/null || echo "0B")
    current_human=$(echo $current_human | sed 's/\([KMGT]B\)B$/\1/')

    if [ $total -eq 0 ]; then
        echo -ne "\r[${filename}] 下载中... | [$(printf "%${PROGRESS_BAR_LENGTH}s" | tr ' ' '-')] 已下载: ${current_human}"
        return
    fi

    local total_human=$(numfmt --to=iec --suffix=B --format="%.1f" $total 2>/dev/null || echo "0B")
    total_human=$(echo $total_human | sed 's/\([KMGT]B\)B$/\1/')
    
    local percentage=$(( current * 100 / total ))
    local filled_length=$(( percentage * PROGRESS_BAR_LENGTH / 100 ))
    local empty_length=$(( PROGRESS_BAR_LENGTH - filled_length ))

    local filled=$(printf "%${filled_length}s" | tr ' ' '#')
    local empty=$(printf "%${empty_length}s" | tr ' ' '-')

    echo -ne "\r[${filename}] ${current_human}/${total_human} | [${filled}${empty}] ${percentage}%"
}

# ===================== 全局进度函数 =====================
show_global_progress() {
    local global_percentage=0
    if [ $TOTAL_FILES -gt 0 ]; then
        global_percentage=$(( DOWNLOADED_FILES * 100 / TOTAL_FILES ))
    fi
    local filled_length=$(( global_percentage * GLOBAL_PROGRESS_LENGTH / 100 ))
    local empty_length=$(( GLOBAL_PROGRESS_LENGTH - filled_length ))

    local filled=$(printf "%${filled_length}s" | tr ' ' '#')
    local empty=$(printf "%${empty_length}s" | tr ' ' '-')

    echo -e "\n📊 全局进度: ${DOWNLOADED_FILES}/${TOTAL_FILES} | [${filled}${empty}] ${global_percentage}% (成功: ${DOWNLOADED_FILES}, 失败: ${FAILED_FILES}, 跳过: ${SKIPPED_FOLDERS}, 重试: ${RETRY_COUNT})"
}

# ===================== 单个文件下载函数 =====================
download_single_file() {
    local alien_path=$1
    local local_path=$2
    local exit_code=0
    local final_size=0

    echo "📥 开始下载: $(basename ${local_path})"

    # 后台下载（--force 强制续传）
    alien_cp --verbose --force --noprompt -t ${TIMEOUT} -T ${THREADS} \
        "alien://${alien_path}" \
        "file:${local_path}" >> ${LOG_FILE} 2>&1 &
    local download_pid=$!

    # 实时监控进度
    local current_size=0
    while kill -0 $download_pid 2>/dev/null; do
        if [ -f "${local_path}" ]; then
            current_size=$(stat -c %s "${local_path}" 2>/dev/null || echo 0)
        fi
        show_progress $current_size 0 "$(basename ${local_path})"
        sleep $CHECK_INTERVAL
    done

    # 等待进程结束
    wait $download_pid
    exit_code=$?

    # 延迟判断文件大小
    echo -e "\n⏳ 等待文件写入完成（延迟${DELAY_CHECK}秒）..."
    sleep $DELAY_CHECK

    # 读取最终大小
    if [ -f "${local_path}" ]; then
        final_size=$(stat -c %s "${local_path}" 2>/dev/null || echo 0)
        local final_human=$(numfmt --to=iec --suffix=B --format="%.1f" $final_size 2>/dev/null || echo "0B")
        echo -e "\r[$(basename ${local_path})] 下载结束 | 大小: ${final_human} (退出码: ${exit_code})"
    else
        final_size=0
        echo -e "\r[$(basename ${local_path})] 无文件生成 (退出码: ${exit_code}) ❌"
    fi

    # 返回结果：0=成功，1=失败
    if [ $exit_code -eq 0 ] && [ $final_size -gt 0 ]; then
        return 0
    else
        return 1
    fi
}

# ===================== 处理单个HY编号的函数 =====================
process_hy_number() {
    local HY_NUM=$1
    
    echo -e "\n====================================="
    echo "======== 处理 hy_${HY_NUM} ========"
    echo "====================================="
    echo "======== 处理 hy_${HY_NUM} ========" >> ${LOG_FILE}
    
    # 构建AOD路径
    AOD_PATH="${ALIEN_BASE_PATH}/hy_${HY_NUM}/AOD"
    
    # 检查AOD路径是否存在
    if ! check_alien_path_exists "${AOD_PATH}"; then
        echo "ℹ️  hy_${HY_NUM}的AOD路径不存在，跳过"
        echo "状态: 跳过 | AOD路径不存在" >> ${LOG_FILE}
        return 1
    fi
    
    # 获取子文件夹列表
    SUB_FOLDERS=()
    while IFS= read -r folder; do
        if [ -n "$folder" ]; then
            SUB_FOLDERS+=("$folder")
        fi
    done < <(get_subfolders "${AOD_PATH}")
    
    if [ ${#SUB_FOLDERS[@]} -eq 0 ]; then
        echo "ℹ️  hy_${HY_NUM}下没有找到数字命名的子文件夹，跳过"
        echo "状态: 跳过 | 无有效子文件夹" >> ${LOG_FILE}
        return 1
    fi
    
    echo "📂 hy_${HY_NUM} 找到 ${#SUB_FOLDERS[@]} 个子文件夹: ${SUB_FOLDERS[*]}"
    
    # 更新总文件数
    TOTAL_FILES=$((TOTAL_FILES + ${#SUB_FOLDERS[@]}))
    
    # 遍历子文件夹
    for SUB_FOLDER in "${SUB_FOLDERS[@]}"; do
        ALIEN_FILE_PATH="${AOD_PATH}/${SUB_FOLDER}/AO2D.root"
        LOCAL_FILE_PATH="${LOCAL_SAVE_DIR}/AO2D_${HY_NUM}_${SUB_FOLDER}.root"

        # 1. 检查AO2D.root文件是否存在
        if ! check_alien_file_exists "${ALIEN_FILE_PATH}"; then
            echo "ℹ️  hy_${HY_NUM}_${SUB_FOLDER} AO2D.root不存在，跳过"
            echo "状态: 跳过 | 文件不存在" >> ${LOG_FILE}
            SKIPPED_FOLDERS=$((SKIPPED_FOLDERS + 1))
            show_global_progress
            continue
        fi

        # 2. 本地文件已存在 → 跳过
        if [ -f "${LOCAL_FILE_PATH}" ]; then
            exist_size=$(stat -c %s "${LOCAL_FILE_PATH}" 2>/dev/null || echo 0)
            exist_size_human=$(numfmt --to=iec --suffix=B --format="%.1f" $exist_size 2>/dev/null || echo "0B")
            echo "✅ hy_${HY_NUM}_${SUB_FOLDER} 本地已存在，跳过（大小: ${exist_size_human}）"
            echo "状态: 跳过 | 本地已存在 | 大小: ${exist_size_human}" >> ${LOG_FILE}
            DOWNLOADED_FILES=$((DOWNLOADED_FILES + 1))
            show_global_progress
            continue
        fi

        # 3. 单文件下载+即时重试
        retry=0
        download_success=0
        echo -e "\n📥 开始处理 hy_${HY_NUM}_${SUB_FOLDER}（最多重试${SINGLE_FILE_RETRY}次）"
        
        while [ "${retry}" -lt "${SINGLE_FILE_RETRY}" ]; do
            if download_single_file "${ALIEN_FILE_PATH}" "${LOCAL_FILE_PATH}"; then
                download_success=1
                break
            else
                retry=$((retry + 1))
                RETRY_COUNT=$((RETRY_COUNT + 1))
                echo "❌ hy_${HY_NUM}_${SUB_FOLDER} 下载失败（重试${retry}/${SINGLE_FILE_RETRY}）"
                echo "状态: 重试 | 次数: ${retry}/${SINGLE_FILE_RETRY}" >> ${LOG_FILE}
                
                # 如果重试次数未用完，删除部分下载的文件
                if [ "${retry}" -lt "${SINGLE_FILE_RETRY}" ] && [ -f "${LOCAL_FILE_PATH}" ]; then
                    rm -f "${LOCAL_FILE_PATH}"
                    echo "  删除不完整文件，准备重试"
                fi
                sleep 5
            fi
        done

        # 4. 最终结果判断
        if [ "${download_success}" -eq 1 ]; then
            file_size=$(stat -c %s "${LOCAL_FILE_PATH}" 2>/dev/null || echo 0)
            file_size_human=$(numfmt --to=iec --suffix=B --format="%.1f" $file_size 2>/dev/null || echo "0B")
            echo "✅ hy_${HY_NUM}_${SUB_FOLDER} 最终下载成功（大小: ${file_size_human}）"
            echo "状态: 成功 | 最终大小: ${file_size_human}" >> ${LOG_FILE}
            DOWNLOADED_FILES=$((DOWNLOADED_FILES + 1))
        else
            echo "❌ hy_${HY_NUM}_${SUB_FOLDER} 重试${SINGLE_FILE_RETRY}次后仍失败，标记为最终失败"
            echo "状态: 最终失败 | 重试次数耗尽" >> ${LOG_FILE}
            FAILED_FILES=$((FAILED_FILES + 1))
            echo "${ALIEN_FILE_PATH}" >> "${LOCAL_SAVE_DIR}/failed_files.txt"
        fi

        # 更新全局进度
        show_global_progress
    done
    
    return 0
}

# ===================== 主程序开始 =====================

# 检查ALIEN环境
check_alien_env

# 创建保存目录
mkdir -p ${LOCAL_SAVE_DIR}

# 初始化日志文件
echo "===== 批量下载开始（边检测边下载版）=====" > ${LOG_FILE}
echo "开始时间: $(date)" >> ${LOG_FILE}
echo "配置：单次超时${TIMEOUT}秒，单文件重试${SINGLE_FILE_RETRY}次，线程数${THREADS}" >> ${LOG_FILE}
echo "ALIEN基础路径: ${ALIEN_BASE_PATH}" >> ${LOG_FILE}
echo "本地保存路径: ${LOCAL_SAVE_DIR}" >> ${LOG_FILE}

# 重置失败列表
> "${LOCAL_SAVE_DIR}/failed_files.txt"

# 初始化统计变量
TOTAL_FILES=0
DOWNLOADED_FILES=0
FAILED_FILES=0
SKIPPED_FOLDERS=0
RETRY_COUNT=0

echo -e "\n🚀 开始批量下载（边检测边下载模式）"

# 遍历HY编号，边检测边下载
for HY_NUM in "${HY_NUMBERS[@]}"; do
    process_hy_number "$HY_NUM"
done

# ===================== 汇总结果 =====================
echo -e "\n\n"
echo "===================== 批量下载结束 ====================="
echo "结束时间: $(date)"
echo "日志文件: ${LOG_FILE}"
echo "失败文件列表: ${LOCAL_SAVE_DIR}/failed_files.txt"

# 计算实际本地文件数
ACTUAL_FILES=$(find ${LOCAL_SAVE_DIR} -name "AO2D_*.root" -type f | wc -l)

echo -e "\n📋 详细统计："
echo "   - 总待处理文件数: ${TOTAL_FILES}"
echo "   - 下载成功数（含本地已存在）: ${DOWNLOADED_FILES}"
echo "   - 最终失败数（重试耗尽）: ${FAILED_FILES}"
echo "   - 跳过数（无文件/无文件夹）: ${SKIPPED_FOLDERS}"
echo "   - 总重试次数: ${RETRY_COUNT}"
echo "   - 实际本地文件数: ${ACTUAL_FILES}"

# 验证统计一致性
EXPECTED_TOTAL=$((DOWNLOADED_FILES + FAILED_FILES + SKIPPED_FOLDERS))
if [ "${EXPECTED_TOTAL}" -eq "${TOTAL_FILES}" ]; then
    echo -e "\n✅ 统计一致性验证通过"
else
    echo -e "\n⚠️  统计一致性警告: 预期总数${TOTAL_FILES}，实际统计${EXPECTED_TOTAL}"
fi

# 写入日志汇总
echo -e "\n===================== 批量下载结束 =====================" >> ${LOG_FILE}
echo "结束时间: $(date)" >> ${LOG_FILE}
echo "总文件数: ${TOTAL_FILES}" >> ${LOG_FILE}
echo "成功数: ${DOWNLOADED_FILES}" >> ${LOG_FILE}
echo "最终失败数: ${FAILED_FILES}" >> ${LOG_FILE}
echo "跳过数: ${SKIPPED_FOLDERS}" >> ${LOG_FILE}
echo "总重试次数: ${RETRY_COUNT}" >> ${LOG_FILE}
echo "实际本地文件数: ${ACTUAL_FILES}" >> ${LOG_FILE}

# 如果有失败文件，给出提示
if [ "${FAILED_FILES}" -gt 0 ]; then
    echo -e "\n⚠️  有 ${FAILED_FILES} 个文件下载失败，可以尝试重新运行脚本（会自动跳过已成功文件）"
    echo "   失败文件列表: ${LOCAL_SAVE_DIR}/failed_files.txt"
fi

echo -e "\n🎉 脚本执行完毕！"