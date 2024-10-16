import aiohttp
import reflex as rx
from Bio.Seq import reverse_complement
from Bio import Align
from Bio.Align import substitution_matrices


Flux_system_prompt = """
你是一个基于Flux.1模型的提示词生成机器人。根据用户的需求，自动生成符合Flux.1格式的绘画提示词。虽然你可以参考提供的模板来学习提示词结构和规律，但你必须具备灵活性来应对各种不同需求。最终输出应仅限提示词，无需任何其他解释或信息。你的回答必须全部使用英语进行回复我！

### **提示词生成逻辑**：

1. **需求解析**：从用户的描述中提取关键信息，包括：
   - 角色：外貌、动作、表情等。
   - 场景：环境、光线、天气等。
   - 风格：艺术风格、情感氛围、配色等。
   - 其他元素：特定物品、背景或特效。

2. **提示词结构规律**：
   - **简洁、精确且具象**：提示词需要简单、清晰地描述核心对象，并包含足够细节以引导生成出符合需求的图像。
   - **灵活多样**：参考下列模板和已有示例，但需根据具体需求生成多样化的提示词，避免固定化或过于依赖模板。
   - **符合Flux.1风格的描述**：提示词必须遵循Flux.1的要求，尽量包含艺术风格、视觉效果、情感氛围的描述，使用与Flux.1模型生成相符的关键词和描述模式。

3. **仅供你参考和学习的几种场景提示词**（你需要学习并灵活调整,"[ ]"中内容视用户问题而定）：
   - **角色表情集**：
场景说明：适合动画或漫画创作者为角色设计多样的表情。这些提示词可以生成展示同一角色在不同情绪下的表情集，涵盖快乐、悲伤、愤怒等多种情感。

提示词：An anime [SUBJECT], animated expression reference sheet, character design, reference sheet, turnaround, lofi style, soft colors, gentle natural linework, key art, range of emotions, happy sad mad scared nervous embarrassed confused neutral, hand drawn, award winning anime, fully clothed

[SUBJECT] character, animation expression reference sheet with several good animation expressions featuring the same character in each one, showing different faces from the same person in a grid pattern: happy sad mad scared nervous embarrassed confused neutral, super minimalist cartoon style flat muted kawaii pastel color palette, soft dreamy backgrounds, cute round character designs, minimalist facial features, retro-futuristic elements, kawaii style, space themes, gentle line work, slightly muted tones, simple geometric shapes, subtle gradients, oversized clothing on characters, whimsical, soft puffy art, pastels, watercolor

   - **全角度角色视图**：
场景说明：当需要从现有角色设计中生成不同角度的全身图时，如正面、侧面和背面，适用于角色设计细化或动画建模。

提示词：A character sheet of [SUBJECT] in different poses and angles, including front view, side view, and back view

   - **80 年代复古风格**：
场景说明：适合希望创造 80 年代复古风格照片效果的艺术家或设计师。这些提示词可以生成带有怀旧感的模糊宝丽来风格照片。

提示词：blurry polaroid of [a simple description of the scene], 1980s.

   - **智能手机内部展示**：
场景说明：适合需要展示智能手机等产品设计的科技博客作者或产品设计师。这些提示词帮助生成展示手机外观和屏幕内容的图像。

提示词：a iphone product image showing the iphone standing and inside the screen the image is shown

   - **双重曝光效果**：
场景说明：适合摄影师或视觉艺术家通过双重曝光技术创造深度和情感表达的艺术作品。

提示词：[Abstract style waterfalls, wildlife] inside the silhouette of a [man]’s head that is a double exposure photograph . Non-representational, colors and shapes, expression of feelings, imaginative, highly detailed

   - **高质感电影海报**：
场景说明：适合需要为电影创建引人注目海报的电影宣传或平面设计师。

提示词：A digital illustration of a movie poster titled [‘Sad Sax: Fury Toad’], [Mad Max] parody poster, featuring [a saxophone-playing toad in a post-apocalyptic desert, with a customized car made of musical instruments], in the background, [a wasteland with other musical vehicle chases], movie title in [a gritty, bold font, dusty and intense color palette].

   - **镜面自拍效果**：
场景说明：适合想要捕捉日常生活瞬间的摄影师或社交媒体用户。

提示词：Phone photo: A woman stands in front of a mirror, capturing a selfie. The image quality is grainy, with a slight blur softening the details. The lighting is dim, casting shadows that obscure her features. [The room is cluttered, with clothes strewn across the bed and an unmade blanket. Her expression is casual, full of concentration], while the old iPhone struggles to focus, giving the photo an authentic, unpolished feel. The mirror shows smudges and fingerprints, adding to the raw, everyday atmosphere of the scene.

   - **像素艺术创作**：
场景说明：适合像素艺术爱好者或复古游戏开发者创造或复刻经典像素风格图像。

提示词：[Anything you want] pixel art style, pixels, pixel art

   - **以上部分场景仅供你学习，一定要学会灵活变通，以适应任何绘画需求**：

4. **Flux.1提示词要点总结**：
   - **简洁精准的主体描述**：明确图像中核心对象的身份或场景。
   - **风格和情感氛围的具体描述**：确保提示词包含艺术风格、光线、配色、以及图像的氛围等信息。
   - **动态与细节的补充**：提示词可包括场景中的动作、情绪、或光影效果等重要细节。
   - **其他更多规律请自己寻找**
---

**问答案例1**：
**用户输入**：一个80年代复古风格的照片。
**你的输出**：A blurry polaroid of a 1980s living room, with vintage furniture, soft pastel tones, and a nostalgic, grainy texture,  The sunlight filters through old curtains, casting long, warm shadows on the wooden floor, 1980s,

**问答案例2**：
**用户输入**：一个赛博朋克风格的夜晚城市背景
**你的输出**：A futuristic cityscape at night, in a cyberpunk style, with neon lights reflecting off wet streets, towering skyscrapers, and a glowing, high-tech atmosphere. Dark shadows contrast with vibrant neon signs, creating a dramatic, dystopian mood
"""

SD_system_prompt = """
作为 Stable Diffusion Prompt 提示词专家，您将从关键词中创建提示，通常来自 Danbooru 等数据库。
提示通常描述图像，使用常见词汇，按重要性排列，并用逗号分隔。避免使用"-"或"."，但可以接受空格和自然语言。避免词汇重复。
为了强调关键词，请将其放在括号中以增加其权重。例如，"(flowers)"将'flowers'的权重增加1.1倍，而"(((flowers)))"将其增加1.331倍。使用"(flowers:1.5)"将'flowers'的权重增加1.5倍。只为重要的标签增加权重。
提示包括三个部分：**前缀** （质量标签+风格词+效果器）+ **主题** （图像的主要焦点）+ **场景** （背景、环境）。
前缀影响图像质量。像"masterpiece"、"best quality"、"4k"这样的标签可以提高图像的细节。像"illustration"、"lensflare"这样的风格词定义图像的风格。像"bestlighting"、"lensflare"、"depthoffield"这样的效果器会影响光照和深度。
主题是图像的主要焦点，如角色或场景。对主题进行详细描述可以确保图像丰富而详细。增加主题的权重以增强其清晰度。对于角色，描述面部、头发、身体、服装、姿势等特征。
场景描述环境。没有场景，图像的背景是平淡的，主题显得过大。某些主题本身包含场景（例如建筑物、风景）。像"花草草地"、"阳光"、"河流"这样的环境词可以丰富场景。你的任务是设计图像生成的提示。请按照以下步骤进行操作：
1. 我会发送给您一个图像场景。需要你生成详细的图像描述
2. 图像描述必须是英文，输出为Positive Prompt。
示例：
我发送：二战时期的护士。
您回复只回复：
A WWII-era nurse in a German uniform, holding a wine bottle and stethoscope, sitting at a table in white attire, with a table in the background, masterpiece, best quality, 4k, illustration style, best lighting, depth of field, detailed character, detailed environment.
"""


async def fetch(session, url, params):
    async with session.get(url, params=params) as response:
        return await response.json()


async def fetch_blat(session, url, params):
    async with session.get(url, params=params) as response:
        return await response.json(content_type="text/html")


async def get_genome_seq(genome, chrom, start, end, reverse) -> dict:
    base_url = "https://api.genome.ucsc.edu/getData/sequence"
    input_params = {
        "genome": genome,
        "chrom": chrom,
        "start": start,
        "end": end,
        "revComp": reverse
    }
    async with aiohttp.ClientSession() as session:
        data = await fetch(session, base_url, input_params)
        return data


async def get_blat_data(userSeq, Type, db, output) -> dict:
    base_url = "https://genome.ucsc.edu/cgi-bin/hgBlat"
    input_params = {
        "userSeq": userSeq,
        "type": Type,
        "db": db,
        "output": output
    }
    async with aiohttp.ClientSession() as session:
        data = await fetch_blat(session, base_url, input_params)
        return data


def parse_coord(coord_string):
    chrom, position = coord_string.split(":")
    start, end = position.split("-")
    return chrom, int(start), int(end)


def parse_form_data(formdata):
    newData = {}
    genome_dict = {
        "Human": "hg38",
        "Zebrafish": "danRer11",
        "Mouse": "mm10"
    }
    newData['genome'] = genome_dict[formdata["species"]]
    chrom, start, end = parse_coord(formdata["coord"])
    newData['chrom'] = chrom
    newData['start'] = start-1
    newData['end'] = end
    if formdata["reverse"] == "3'-5'":
        newData['reverse'] = 1
    else:
        newData['reverse'] = 0
    return newData


def parse_blat_form(formdata):
    newData = {}
    genome_dict = {
        "Human": "hg38",
        "Zebrafish": "danRer11",
        "Mouse": "mm10"
    }
    newData['db'] = genome_dict[formdata["genome"]]
    newData['Type'] = formdata['query_type']
    newData['userSeq'] = formdata['seq']
    newData['output'] = "json"
    return newData


def show_blat(blat_result_row: list):
    return rx.table.row(
        rx.table.cell(blat_result_row[0]),
        rx.table.cell(blat_result_row[1]),
        rx.table.cell(blat_result_row[2]),
        rx.table.cell(blat_result_row[3]),
        rx.table.cell(blat_result_row[4]),
        rx.table.cell(blat_result_row[5]),
        rx.table.cell(blat_result_row[6]),
        rx.table.cell(blat_result_row[7])
    )


def rev_comp(seq):
    return reverse_complement(seq)


def to_uppercase(seq):
    return seq.upper()


def to_lowercase(seq):
    return seq.lower()


def global_alignment(ref, query, gap_open, gap_extend, end_gap_open, end_gap_extend, end_gap):
    matrix = substitution_matrices.load("NUC.4.4")
    if not end_gap:
        aligner = Align.PairwiseAligner(
            substitution_matrix=matrix,
            mode="global",
            internal_open_gap_score=0-gap_open,
            internal_extend_gap_score=0-gap_extend,
            end_gap_score=0
        )
    else:
        aligner = Align.PairwiseAligner(
            substitution_matrix=matrix,
            mode="global",
            internal_open_gap_score=0-gap_open,
            internal_extend_gap_score=0-gap_extend,
            end_open_gap_score=0-end_gap_open,
            end_extend_gap_score=0-end_gap_extend
        )
    ref = ref.strip().upper()
    query = query.strip().upper()
    alignments = aligner.align(ref, query)
    return alignments[0].format(), alignments[0].score


def render_codes(code_row: str) -> rx.Component:
    return rx.text(code_row, class_name="font-mono")

#################################################################################################################
############################ AI ART PAGE##########################################################################
#################################################################################################################





def get_prompt(client, prompt, model="gpt-4o", modeltype="flux"):
    if modeltype == "flux":
        completion = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": Flux_system_prompt},
                {
                    "role": "user",
                    "content": prompt
                }
            ]
        )
        return completion.choices[0].message.content
    else:
        completion = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": SD_system_prompt},
                {
                    "role": "user",
                    "content": prompt
                }
            ]
        )
        return completion.choices[0].message.content


async def get_picture_SC(client, model, api_token, prompt, negtive_prompt, image_size, steps, scale, translate, chatmodel="gpt-4o"):
    base_url = "https://api.siliconflow.cn/v1/image/generations"
    headers = {
        "Authorization": f"Bearer {api_token}",
        "Content-Type": "application/json"
    }
    model_dict = {
        "FLUX.1-schnell": "black-forest-labs/FLUX.1-schnell",
        "FLUX.1-dev": "black-forest-labs/FLUX.1-dev",
        "stable-diffusion-3-medium": "stabilityai/stable-diffusion-3-medium",
        "stable-diffusion-xl-base-1.0": "stabilityai/stable-diffusion-xl-base-1.0"
    }
    if translate:
        if model == "FLUX.1-schnell" or model == "FLUX.1-dev":
            NewPrompts = get_prompt(client,prompt,model=chatmodel,modeltype="flux")
        else:
            NewPrompts = get_prompt(client, prompt, model=chatmodel,modeltype="sd")
    else:
        NewPrompts = prompt
    if model == "FLUX.1-schnell":
        payload = {
            "model": model_dict[model],
            "prompt": NewPrompts,
            "image_size": image_size,
        }
    elif model == "FLUX.1-dev":
        payload = {
            "model": model_dict[model],
            "prompt": NewPrompts,
            "image_size": image_size,
            "num_inference_steps": steps,
        }
    else:
        payload = {
            "model": model_dict[model],
            "prompt": NewPrompts,
            "negative_prompt": negtive_prompt,
            "image_size": image_size,
            "batch_size": 1,
            "num_inference_steps": steps,
            "guidance_scale": scale
        }
    async with aiohttp.ClientSession() as session:
        async with session.post(base_url, json=payload, headers=headers) as response:
            if response.status == 200:
                data = await response.json()
                return data['images'][0]['url'], NewPrompts
            else:
                return None, NewPrompts

#################################################################################################################
