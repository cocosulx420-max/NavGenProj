--!strict
-- NavGen.Boundary — clean outlines from the local grids, by DESIGN.md's method.
--
-- The input is LocalGrid's per-part, part-aligned cell masks. The method is
-- DESIGN.md steps 3-5, 7 and 8: chain the boundary edges of a mask, grow a
-- best-fit line along the boundary cells while the MAXIMUM residual stays under
-- a cell, and take corners as the places where that fit failed.
--
-- THIS MODULE NEVER ASKS A PART ANYTHING. No CFrames, no face planes, no sizes,
-- no raycasts, no overlap queries. It reads cell centres and nothing else. That
-- is DESIGN.md's central rule and it is not an aesthetic one: a Union or a
-- MeshPart has no readable planar face, and interpenetrating parts report faces
-- that are not surfaces. An earlier version of this module took its lines from
-- the side planes of whichever part killed each cell. It was more accurate on
-- SmallMap -- and it is exactly the approach DESIGN.md records as having given
-- this project hell on real maps, so it is gone.
--
-- WHY THE FIT IS CHEAP HERE, WHICH IS THE POINT OF THE LOCAL GRIDS. On a
-- world-aligned raster every rotated part staircases, so the fit is doing heavy
-- reconstruction everywhere and its tolerance is load-bearing everywhere. On a
-- part-aligned grid the host's own rim already lies along whole lattice lines,
-- so the fit reproduces it with zero residual and has nothing to undo. What is
-- left for the fit is the genuinely hard case: the footprint of some OTHER part
-- crossing this one's lattice at an angle. Same method as DESIGN.md, applied to
-- a fraction of the edges, which is why the tolerances stop fighting.
--
-- SAFETY. The only morphology is erosion. Lines are biased inward (step 5) and
-- polygons are fitted to cell centres, so the boundary sits at or inside the
-- walkable cells and never outside them. Erosion can remove walkable ground but
-- can never invent connectivity, so no amount of it welds two rooms through a
-- wall. Nothing here dilates.

local Boundary = {}

export type Config = {
	-- Fallback grids only. A world-aligned grid can hold two surfaces at once
	-- (a mesh with a ledge), so its trace needs to know what counts as a cliff.
	-- A block grid is a single plane and never consults this.
	stepTol: number?,

	-- STEP 4. The maximum perpendicular distance a boundary cell centre may sit
	-- from its segment's line. One cell: anything the mask can express as
	-- "straight" is straight.
	fitTol: number?,

	-- STEP 8. Runs shorter than this are absorbed, and runs whose directions
	-- agree to within collinearDeg are merged -- both BEFORE corners are
	-- intersected, because short runs are what make an intersection unstable.
	minSegLen: number?,
	collinearDeg: number?,

	-- STEP 7. An acute corner throws two lines' intersection arbitrarily far
	-- out -- the classic miter spike. Past this distance from the corner it
	-- replaces, bevel across instead.
	miterLimit: number?,
}

local DEFAULT = {
	stepTol = 2.2,
	fitTol = 1.0,
	minSegLen = 1.0,
	collinearDeg = 5,
	miterLimit = 3.0,
}

local function merged(cfg): any
	local c = {}
	for k, v in pairs(DEFAULT) do c[k] = v end
	if cfg then for k, v in pairs(cfg) do if v ~= nil then c[k] = v end end end
	return c
end

--------------------------------------------------------------------------
-- 2D helpers, in the grid's own face coordinates
--------------------------------------------------------------------------

type P2 = { x: number, z: number }

local function sub(a: P2, b: P2): P2 return { x = a.x - b.x, z = a.z - b.z } end
local function dot(a: P2, b: P2): number return a.x * b.x + a.z * b.z end
local function len(a: P2): number return math.sqrt(a.x * a.x + a.z * a.z) end

-- Interior lies to the LEFT of travel, so this is the outward side.
local function outwardOf(d: P2): P2
	local m = len(d)
	if m < 1e-9 then return { x = 0, z = 0 } end
	return { x = d.z / m, z = -d.x / m }
end

-- Total least squares: the principal axis of the point set. NOT ordinary least
-- squares -- a boundary run can be near-vertical in the grid's frame, where
-- fitting one coordinate as a function of the other blows up.
--
-- Takes an explicit list of point INDICES rather than a range: the run that
-- straddles the loop's arbitrary start point wraps, and a range cannot express
-- that.
local function fitLine(pts: {P2}, idx: {number}): (P2, P2)
	local n = #idx
	local sx, sz = 0, 0
	for _, i in ipairs(idx) do sx += pts[i].x; sz += pts[i].z end
	local cx, cz = sx / n, sz / n
	local sxx, szz, sxz = 0, 0, 0
	for _, i in ipairs(idx) do
		local dx, dz = pts[i].x - cx, pts[i].z - cz
		sxx += dx * dx; szz += dz * dz; sxz += dx * dz
	end
	-- larger eigenvalue of the 2x2 covariance, closed form
	local tr, det = sxx + szz, sxx * szz - sxz * sxz
	local disc = math.max(tr * tr * 0.25 - det, 0)
	local lam = tr * 0.5 + math.sqrt(disc)
	local dx, dz
	if math.abs(sxz) > 1e-12 then
		dx, dz = lam - szz, sxz
	elseif sxx >= szz then
		dx, dz = 1, 0
	else
		dx, dz = 0, 1
	end
	-- ORIENT ALONG TRAVEL. A principal axis has no sign -- PCA is just as happy
	-- to hand back the reverse direction -- but everything downstream reads the
	-- sign as meaning something. outwardOf() assumes interior-on-the-left, which
	-- is only true if dir runs the way the loop was walked; a flipped segment
	-- gets its "outward" normal pointing INTO the floor, and the inward bias
	-- then pushes the line the wrong way, straight through the wall.
	local sxs = pts[idx[n]].x - pts[idx[1]].x
	local szs = pts[idx[n]].z - pts[idx[1]].z
	if dx * sxs + dz * szs < 0 then dx, dz = -dx, -dz end
	local m = math.sqrt(dx * dx + dz * dz)
	if m < 1e-12 then dx, dz, m = 1, 0, 1 end
	return { x = cx, z = cz }, { x = dx / m, z = dz / m }
end

local function maxResidual(pts: {P2}, idx: {number}, c: P2, d: P2): number
	local nx, nz = -d.z, d.x
	local worst = 0
	for _, i in ipairs(idx) do
		local r = math.abs((pts[i].x - c.x) * nx + (pts[i].z - c.z) * nz)
		if r > worst then worst = r end
	end
	return worst
end

--------------------------------------------------------------------------
-- STEP 3 — trace the contour of a cell mask
--
-- Boundary EDGES are chained rather than boundary cells walked: the boundary of
-- a set of cells is a closed loop by construction, so there is no open end to
-- chase and no tolerance involved. The ordered cells the fit needs fall out of
-- the edge order.
--
-- `cliff(a, b)` reports that two lattice-adjacent cells are not continuous
-- ground even though both are in the mask. On a block grid it is never true;
-- on a fallback grid it is what stops a polygon spanning a mesh's ledge.
--------------------------------------------------------------------------

local DIRS = { { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } }

-- lattice corners of cell (u,v)'s edge in direction di, oriented so the loop
-- runs with the interior on its left
local function cornersOf(u: number, v: number, di: number)
	if di == 1 then return { u + 1, v }, { u + 1, v + 1 }
	elseif di == 2 then return { u + 1, v + 1 }, { u, v + 1 }
	elseif di == 3 then return { u, v + 1 }, { u, v }
	else return { u, v }, { u + 1, v } end
end

local function traceMask(cells: {any}, index: {[string]: any}, cliff: ((any, any) -> boolean)?)
	local segs: {any} = {}
	local byStart: {[string]: {any}} = {}
	for _, cell in ipairs(cells) do
		for di, d in ipairs(DIRS) do
			local nb = index[(cell.ui + d[1]) .. ":" .. (cell.vi + d[2])]
			if not nb or (cliff and cliff(cell, nb)) then
				local a, b = cornersOf(cell.ui, cell.vi, di)
				local s = { a = a, b = b, cell = cell }
				segs[#segs + 1] = s
				local k = a[1] .. "," .. a[2]
				local bucket = byStart[k]
				if not bucket then bucket = {}; byStart[k] = bucket end
				bucket[#bucket + 1] = s
			end
		end
	end

	-- CHAINING AT A JUNCTION is not a free choice. A cliff between two cells
	-- that are both in the mask emits an edge from EACH side, on the same
	-- lattice edge in opposite directions, because the upper rim and the lower
	-- rim are two different boundaries that coincide in plan. The outline
	-- therefore has zero-width slits whose ends are vertices with four edges
	-- leaving them, and taking whichever candidate came first there splices one
	-- rim onto the other and drops degenerate slivers out of the walk.
	--
	-- Standard face traversal of an embedded planar graph fixes it: arriving
	-- along an edge, leave on the next edge CLOCKWISE from it. With the interior
	-- kept on the left, that walks each rim whole and turns a slit around at its
	-- tip.
	local function angle(s: any): number
		return math.atan2(s.b[2] - s.a[2], s.b[1] - s.a[1])
	end
	local TAU = math.pi * 2
	local used: {[any]: boolean} = {}
	local loops: {{any}} = {}
	for _, s0 in ipairs(segs) do
		if not used[s0] then
			local loop, cur = {}, s0
			while cur and not used[cur] do
				used[cur] = true
				loop[#loop + 1] = cur
				local cands = byStart[cur.b[1] .. "," .. cur.b[2]]
				local nxt, bestTurn = nil, math.huge
				if cands then
					local back = angle(cur) + math.pi
					for _, s in ipairs(cands) do
						if not used[s] then
							local turn = (back - angle(s)) % TAU
							-- straight back is the full turn, taken only when
							-- nothing else is left: the tip of a slit
							if turn <= 1e-9 then turn = TAU end
							if turn < bestTurn then nxt, bestTurn = s, turn end
						end
					end
				end
				cur = nxt
			end
			if #loop >= 4 then loops[#loops + 1] = loop end
		end
	end
	return loops
end

--------------------------------------------------------------------------
-- STEP 4 — greedy line fit. THIS IS WHERE THE STAIRCASE DIES.
--
-- Walk the loop maintaining a best-fit line through the cells accepted so far.
-- While the MAXIMUM perpendicular residual stays under a cell, keep extending.
-- When it exceeds, close the segment and start fresh from that cell.
--
-- The cells of a foreign part's footprint crossing this grid at an angle all
-- sit within a cell of one straight line at the true angle, so they collapse
-- into a single segment and the fit recovers the angle without being told it.
--
-- CORNERS ARE WHERE THE FIT FAILS. They are a byproduct, not a prerequisite --
-- which is the piece that killed every attempt that tried to detect them
-- against real geometry instead.
--------------------------------------------------------------------------

local function segmentLoop(pts: {P2}, c: any)
	local n = #pts
	local segs: {any} = {}
	local cur = { 1 }
	for i = 2, n do
		local trial = table.clone(cur)
		trial[#trial + 1] = i
		local cen, dir = fitLine(pts, trial)
		-- MAXIMUM residual, never the average: an average lets a shallow corner
		-- hide inside a long run, which is precisely the corner worth keeping.
		if maxResidual(pts, trial, cen, dir) > c.fitTol then
			-- Close the run and start the next one AT this point, so the two
			-- segments share it.
			local ce, de = fitLine(pts, cur)
			segs[#segs + 1] = { idx = cur, cen = ce, dir = de }
			cur = { i - 1, i }
		else
			cur = trial
		end
	end
	local ce, de = fitLine(pts, cur)
	segs[#segs + 1] = { idx = cur, cen = ce, dir = de }
	return segs
end

local function spanLength(pts: {P2}, idx: {number}): number
	return len(sub(pts[idx[#idx]], pts[idx[1]]))
end

--------------------------------------------------------------------------
-- STEP 8 — clean up, and all of it BEFORE any corner is intersected
--------------------------------------------------------------------------

local function mergeSegments(pts: {P2}, segs: {any}, c: any)
	local cosLim = math.cos(math.rad(c.collinearDeg))

	local function refit(a: any, b: any): any?
		local idx = table.clone(a.idx)
		local seen: {[number]: boolean} = {}
		for _, i in ipairs(idx) do seen[i] = true end
		for _, i in ipairs(b.idx) do
			if not seen[i] then idx[#idx + 1] = i; seen[i] = true end
		end
		local cen, dir = fitLine(pts, idx)
		if maxResidual(pts, idx, cen, dir) > c.fitTol then return nil end
		return { idx = idx, cen = cen, dir = dir }
	end

	-- SEAM. The walk's start point on a closed loop is arbitrary, so a straight
	-- run that happens to straddle it is always split in two -- a spurious
	-- corner at exactly the place where nothing happened.
	if #segs >= 2 then
		local first, last = segs[1], segs[#segs]
		if dot(first.dir, last.dir) >= cosLim then
			local m = refit(last, first)
			if m then
				segs[1] = m
				segs[#segs] = nil
			end
		end
	end

	-- Near-collinear and too-short runs, absorbed until nothing changes.
	--
	-- This is also the rule that fixes a staircase built out of SEPARATE PARTS.
	-- A flight here is nine stacked blocks whose ends are flush, so the foot of
	-- the flight is one straight edge crossed by eight part seams. Each seam
	-- puts a one-cell jog in the mask, the fit closes a run at each, and without
	-- this merge the result is eight spurious corners along a straight line.
	local changed = true
	while changed and #segs > 3 do
		changed = false
		for k = 1, #segs do
			local a, b = segs[k], segs[(k % #segs) + 1]
			if a == b then break end
			-- SIGNED, not abs. Directions are oriented along the walk, so an abs
			-- test calls a 180-degree reversal "collinear" -- and the two long
			-- sides of a thin wall's ring are exactly that, sitting within fitTol
			-- of each other because the wall is thinner than the tolerance. They
			-- would merge and the ring would collapse.
			local nearly = dot(a.dir, b.dir) >= cosLim
			local tiny = spanLength(pts, a.idx) < c.minSegLen
				or spanLength(pts, b.idx) < c.minSegLen
			if nearly or tiny then
				local m = refit(a, b)
				if m then
					segs[k] = m
					table.remove(segs, (k % #segs) + 1)
					changed = true
					break
				end
			end
		end
	end
	return segs
end

local function signed2(pts: {P2}): number
	local a2 = 0
	for i = 1, #pts do
		local p, q = pts[i], pts[i % #pts + 1]
		a2 += p.x * q.z - q.x * p.z
	end
	return a2 * 0.5
end

--------------------------------------------------------------------------
-- One grid -> its rings, in world space
--------------------------------------------------------------------------

local function ringsOfGrid(g: any, c: any, stats: any)
	local isBlock = not g.fallback and g.n ~= nil
	local step = g.step

	local toWorld
	if isBlock then
		-- A point in the host's face coordinates is a point ON the host's face
		-- plane, so heights are exact by construction: nothing is sampled from a
		-- neighbouring cell and a vertex cannot pick up a height from the wrong
		-- side of a cliff.
		toWorld = function(p: P2): Vector3
			return g.origin + g.u * p.x + g.v * p.z
		end
	else
		-- Fallback grid: coordinates are world XZ and the surface is not one
		-- plane, so a point takes the height of the highest cell touching it.
		toWorld = function(p: P2): Vector3
			local iu, iv = math.floor(p.x / step), math.floor(p.z / step)
			local bestY = nil
			for du = -1, 1 do
				for dv = -1, 1 do
					local cell = g.index[(iu + du) .. ":" .. (iv + dv)]
					if cell and (not bestY or cell.pos.Y > bestY) then bestY = cell.pos.Y end
				end
			end
			return Vector3.new(p.x, bestY or 0, p.z)
		end
	end

	local cliff = nil
	if not isBlock then
		cliff = function(a: any, b: any): boolean
			return math.abs(a.pos.Y - b.pos.Y) > c.stepTol
		end
	end

	-- A ring is built in the grid's 2D face coordinates and only converted to
	-- world once outer-vs-hole is known, because a ring that could not be fitted
	-- is treated differently depending on which it is.
	local rings: {any} = {}
	for _, loop in ipairs(traceMask(g.cells, g.index, cliff)) do
		-- Ordered boundary cell CENTRES. A corner cell contributes two edges and
		-- the duplicate carries no information, so collapse it.
		local pts: {P2} = {}
		local lastCell = nil
		for _, s in ipairs(loop) do
			if s.cell ~= lastCell then
				pts[#pts + 1] = { x = (s.cell.ui + 0.5) * step, z = (s.cell.vi + 0.5) * step }
				lastCell = s.cell
			end
		end
		if #pts < 3 then continue end

		-- A RING TOO THIN TO SURVIVE THE FIT. DESIGN.md step 8 warns about this
		-- and the old implementation hit it on 165 of 169 holes: the two long
		-- sides of a strip narrower than fitTol sit within the tolerance of each
		-- other, so the fit walks straight round the end without ever failing and
		-- the whole ring collapses to one or two segments.
		--
		-- It is not a corner case here. This map's stair steps are 35x2 cell
		-- strips -- 36 of them, 3.6% of every walkable cell on the map -- and
		-- silently dropping them is the one operation in this module that can
		-- make real ground disappear. So a ring that cannot be fitted keeps its
		-- raw lattice outline instead. Jagged, but present.
		local segs = mergeSegments(pts, segmentLoop(pts, c), c)
		stats.rawSegments += #segs
		if #segs < 3 then
			stats.rawRings += 1
			rings[#rings + 1] = { pts2 = pts, area = signed2(pts), raw = true }
			continue
		end

		----------------------------------------------------------------
		-- STEP 5 — bias the fit inward
		--
		-- A fit can sit outward of the cells it was fitted to, which hands back
		-- ground that is not walkable. Translate each line inward until no
		-- accepted cell centre lies outward of it. Without this the safety
		-- guarantee is probabilistic; with it, it is exact.
		----------------------------------------------------------------
		local lines: {any} = {}
		for _, sg in ipairs(segs) do
			local nrm = outwardOf(sg.dir)
			local cval = -math.huge
			for _, i in ipairs(sg.idx) do
				local d = dot(pts[i], nrm)
				if d > cval then cval = d end
			end
			lines[#lines + 1] = { n = nrm, c = cval, anchor = pts[sg.idx[#sg.idx]] }
		end

		----------------------------------------------------------------
		-- STEP 7 — corners are the intersections of adjacent lines
		----------------------------------------------------------------
		local verts: {P2} = {}
		local nL = #lines
		for i = 1, nL do
			local l1, l2 = lines[i], lines[(i % nL) + 1]
			local det = l1.n.x * l2.n.z - l1.n.z * l2.n.x
			local anchor = l1.anchor
			if math.abs(det) < 1e-6 then
				-- near-parallel: fall back to the foot of the anchor on l1
				stats.unstableCorners += 1
				local s = l1.c - dot(anchor, l1.n)
				verts[#verts + 1] = { x = anchor.x + l1.n.x * s, z = anchor.z + l1.n.z * s }
			else
				local px = (l1.c * l2.n.z - l2.c * l1.n.z) / det
				local pz = (l1.n.x * l2.c - l2.n.x * l1.c) / det
				if len(sub({ x = px, z = pz }, anchor)) > c.miterLimit then
					-- MITER LIMIT. An acute corner throws the intersection
					-- arbitrarily far out; bevel across it instead.
					stats.bevels += 1
					local s1 = l1.c - dot(anchor, l1.n)
					verts[#verts + 1] = { x = anchor.x + l1.n.x * s1, z = anchor.z + l1.n.z * s1 }
					local s2 = l2.c - dot(anchor, l2.n)
					verts[#verts + 1] = { x = anchor.x + l2.n.x * s2, z = anchor.z + l2.n.z * s2 }
				else
					verts[#verts + 1] = { x = px, z = pz }
				end
			end
		end
		if #verts < 3 then
			stats.rawRings += 1
			rings[#rings + 1] = { pts2 = pts, area = signed2(pts), raw = true }
			continue
		end
		rings[#rings + 1] = { pts2 = verts, area = signed2(verts) }
	end

	-- Outer is the largest by magnitude. Magnitude, never the sign: this project
	-- has been bitten by a handedness assumption before.
	local outer = rings[1]
	for _, r in ipairs(rings) do
		if math.abs(r.area) > math.abs(outer.area) then outer = r end
	end

	for _, r in ipairs(rings) do
		-- An unfitted OUTER ring stands as it is: cell centres lie inside the
		-- walkable cells, so the polygon is conservative already. An unfitted
		-- HOLE is the opposite -- cell centres sit half a cell INTO the obstacle
		-- the hole represents, so the hole would come out too small and hand
		-- back ground that is not there. Push it out by half a cell. An obstacle
		-- that is slightly too big is safe; one that is too small is not.
		if r.raw and r ~= outer then
			local cx, cz = 0, 0
			for _, p in ipairs(r.pts2) do cx += p.x; cz += p.z end
			cx, cz = cx / #r.pts2, cz / #r.pts2
			for _, p in ipairs(r.pts2) do
				local dx, dz = p.x - cx, p.z - cz
				local m = math.sqrt(dx * dx + dz * dz)
				if m > 1e-6 then
					p.x += dx / m * step * 0.5
					p.z += dz / m * step * 0.5
				end
			end
		end
		local world = table.create(#r.pts2)
		for i, p in ipairs(r.pts2) do world[i] = toWorld(p) end
		r.verts = world
		r.outer = (r == outer)
		stats.segments += #world
	end
	return rings
end

--------------------------------------------------------------------------
-- Boundary.fromLocal — the entry point
--------------------------------------------------------------------------

function Boundary.fromLocal(localData: any, cfg: Config?)
	local c = merged(cfg)
	local t0 = os.clock()

	local regions: {any} = {}
	local stats = {
		parts = 0, block = 0, fallback = 0,
		regions = 0, holes = 0, verts = 0, emptyGrids = 0,
		rawSegments = 0, segments = 0, bevels = 0, unstableCorners = 0,
		-- rings the fit could not resolve, kept as their raw lattice outline
		rawRings = 0,
	}

	for part, g in pairs(localData.grids) do
		stats.parts += 1
		if g.fallback then stats.fallback += 1 else stats.block += 1 end
		local rings = ringsOfGrid(g, c, stats)
		if #rings == 0 then
			stats.emptyGrids += 1
			continue
		end
		local outer = rings[1]
		for _, r in ipairs(rings) do
			if r.outer then outer = r end
		end
		local region = {
			part = part,
			fallback = g.fallback,
			verts = outer.verts,
			area = math.abs(outer.area),
			cells = #g.cells,
			dead = #g.dead,
			holes = {},
		}
		for _, r in ipairs(rings) do
			if r ~= outer then
				region.holes[#region.holes + 1] = { verts = r.verts, area = math.abs(r.area) }
			end
		end
		stats.verts += #region.verts
		for _, h in ipairs(region.holes) do stats.verts += #h.verts end
		stats.regions += 1
		stats.holes += #region.holes
		regions[#regions + 1] = region
	end

	stats.seconds = os.clock() - t0
	return { regions = regions, stats = stats, config = c }
end

--------------------------------------------------------------------------

function Boundary.visualize(res: any, parent: Instance?)
	local root = parent or workspace
	local old = root:FindFirstChild("NavGen_Boundary")
	if old then old:Destroy() end
	local folder = Instance.new("Folder")
	folder.Name = "NavGen_Boundary"
	folder.Parent = root

	local function bar(a: Vector3, b: Vector3, colour: Color3, thick: number, into: Instance)
		local d = b - a
		if d.Magnitude < 1e-4 then return end
		local p = Instance.new("Part")
		p.Anchored = true; p.CanCollide = false; p.CanQuery = false; p.CanTouch = false
		p.Material = Enum.Material.Neon
		p.Color = colour
		p.Size = Vector3.new(thick, thick, d.Magnitude)
		p.CFrame = CFrame.lookAt(a + d * 0.5, b) + Vector3.new(0, 0.08, 0)
		p.Parent = into
	end

	-- Outer rings red, holes cyan, and a FALLBACK part's outline in amber --
	-- those are the world-aligned ones, where the fit is doing the heavy
	-- reconstruction, so they should be identifiable without counting anything.
	for ri, r in ipairs(res.regions) do
		local sub = Instance.new("Folder")
		sub.Name = string.format("R%d_%s_c%d_h%d%s", ri, r.part.Name, r.cells, #r.holes,
			r.fallback and "_FALLBACK" or "")
		sub.Parent = folder
		local outerCol = r.fallback and Color3.fromRGB(255, 170, 60) or Color3.fromRGB(255, 80, 80)
		local v = r.verts
		for i = 1, #v do
			bar(v[i], v[i % #v + 1], outerCol, 0.18, sub)
		end
		for _, hr in ipairs(r.holes) do
			local hv = hr.verts
			for i = 1, #hv do
				bar(hv[i], hv[i % #hv + 1], Color3.fromRGB(80, 200, 255), 0.18, sub)
			end
		end
	end
	return folder
end

return Boundary
